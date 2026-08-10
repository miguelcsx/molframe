//! Leakage-aware manifest splitting.

use std::collections::BTreeMap;

use crate::numeric::{f64_to_usize, u64_to_usize, usize_to_f64, usize_to_u64};

use pdbiox_seq::{Scoring, global};

use super::{Dataset, DatasetError};

const RATIO_TOLERANCE: f64 = 1.0e-9;

/// Requested train, validation and test proportions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitRatios {
    /// Training proportion.
    pub train: f64,
    /// Validation proportion.
    pub validation: f64,
    /// Test proportion.
    pub test: f64,
}

impl SplitRatios {
    /// Creates finite, non-negative ratios that sum to one.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError::InvalidRatios`] for an invalid distribution.
    pub fn new(train: f64, validation: f64, test: f64) -> Result<Self, DatasetError> {
        let values = [train, validation, test];
        if values
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
            || (values.iter().sum::<f64>() - 1.0).abs() > RATIO_TOLERANCE
            || train == 0.0
        {
            return Err(DatasetError::InvalidRatios);
        }
        Ok(Self {
            train,
            validation,
            test,
        })
    }

    const fn values(self) -> [f64; 3] {
        [self.train, self.validation, self.test]
    }
}

/// Manifest-only split strategy.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum SplitStrategy {
    /// Single-linkage groups from globally aligned sequence identity.
    SequenceIdentity {
        /// Minimum identity connecting two entries, in `[0, 1]`.
        threshold: f64,
    },
    /// Keep every precomputed structural cluster in one partition.
    StructuralCluster,
    /// Earlier deposition dates train; later dates validate and test.
    Temporal,
    /// Entry-level deterministic shuffle. This does not prevent leakage.
    Random {
        /// Explicit reproducibility seed.
        seed: u64,
    },
}

/// Complete split request.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitOptions {
    /// Grouping or ordering policy.
    pub strategy: SplitStrategy,
    /// Desired proportions.
    pub ratios: SplitRatios,
}

/// Non-fatal scientific warning attached to a split.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum DatasetWarning {
    /// Random splitting makes no redundancy guarantee.
    RandomSplitMayLeak,
}

/// Three immutable dataset views and warnings produced by the policy.
#[derive(Clone, Debug, PartialEq)]
pub struct DatasetSplit {
    /// Training subset.
    pub train: Dataset,
    /// Validation subset.
    pub validation: Dataset,
    /// Test subset.
    pub test: Dataset,
    /// Explicit policy warnings.
    pub warnings: Vec<DatasetWarning>,
}

pub(crate) fn split(
    dataset: &Dataset,
    options: &SplitOptions,
) -> Result<DatasetSplit, DatasetError> {
    let (partitions, warnings) = match options.strategy {
        SplitStrategy::SequenceIdentity { threshold } => {
            if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
                return Err(DatasetError::InvalidThreshold);
            }
            (
                assign_groups(sequence_groups(dataset, threshold)?, options.ratios),
                Vec::new(),
            )
        }
        SplitStrategy::StructuralCluster => (
            assign_groups(structural_groups(dataset)?, options.ratios),
            Vec::new(),
        ),
        SplitStrategy::Temporal => (temporal(dataset, options.ratios)?, Vec::new()),
        SplitStrategy::Random { seed } => (
            random(dataset, options.ratios, seed),
            vec![DatasetWarning::RandomSplitMayLeak],
        ),
    };
    Ok(DatasetSplit {
        train: dataset.subset(partitions[0].clone()),
        validation: dataset.subset(partitions[1].clone()),
        test: dataset.subset(partitions[2].clone()),
        warnings,
    })
}

fn structural_groups(dataset: &Dataset) -> Result<Vec<Vec<usize>>, DatasetError> {
    let mut groups: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for index in dataset.absolute_indices() {
        let entry = dataset.entry(*index);
        let cluster = entry.structure_cluster.as_deref().ok_or_else(|| {
            DatasetError::MissingSplitMetadata {
                id: entry.id.clone(),
                field: "structure_cluster",
            }
        })?;
        groups.entry(cluster).or_default().push(*index);
    }
    Ok(groups.into_values().collect())
}

fn sequence_groups(dataset: &Dataset, threshold: f64) -> Result<Vec<Vec<usize>>, DatasetError> {
    let indices = dataset.absolute_indices();
    let sequences = indices
        .iter()
        .map(|index| {
            let entry = dataset.entry(*index);
            entry
                .sequence
                .as_deref()
                .ok_or_else(|| DatasetError::MissingSplitMetadata {
                    id: entry.id.clone(),
                    field: "sequence",
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut components = DisjointSet::new(indices.len());
    for left in 0..sequences.len() {
        for right in left + 1..sequences.len() {
            if sequence_identity(sequences[left].as_bytes(), sequences[right].as_bytes())?
                >= threshold
            {
                components.union(left, right);
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (local, absolute) in indices.iter().enumerate() {
        groups
            .entry(components.find(local))
            .or_default()
            .push(*absolute);
    }
    Ok(groups.into_values().collect())
}

fn sequence_identity(left: &[u8], right: &[u8]) -> Result<f64, DatasetError> {
    let alignment = global(left, right, Scoring::simple())?;
    let aligned = alignment
        .columns
        .iter()
        .filter(|column| column.left.is_some() && column.right.is_some())
        .count();
    if aligned == 0 {
        return Ok(0.0);
    }
    let identical = alignment
        .columns
        .iter()
        .filter(|column| {
            column
                .left
                .zip(column.right)
                .is_some_and(|(a, b)| left[a].eq_ignore_ascii_case(&right[b]))
        })
        .count();
    Ok(usize_to_f64(identical) / usize_to_f64(aligned))
}

fn temporal(dataset: &Dataset, ratios: SplitRatios) -> Result<[Vec<usize>; 3], DatasetError> {
    let mut dated = dataset
        .absolute_indices()
        .iter()
        .map(|index| {
            let entry = dataset.entry(*index);
            entry
                .deposition_date
                .as_deref()
                .map(|date| (date, *index))
                .ok_or_else(|| DatasetError::MissingSplitMetadata {
                    id: entry.id.clone(),
                    field: "deposition_date",
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    dated.sort_unstable();
    let ordered = dated
        .into_iter()
        .map(|(_, index)| index)
        .collect::<Vec<_>>();
    Ok(slice_by_ratios(&ordered, ratios))
}

fn random(dataset: &Dataset, ratios: SplitRatios, seed: u64) -> [Vec<usize>; 3] {
    let mut indices = dataset.absolute_indices().to_vec();
    let mut random = SplitMix64(seed);
    for end in (1..indices.len()).rev() {
        let selected = u64_to_usize(random.next() % (usize_to_u64(end) + 1));
        indices.swap(end, selected);
    }
    slice_by_ratios(&indices, ratios)
}

fn slice_by_ratios(indices: &[usize], ratios: SplitRatios) -> [Vec<usize>; 3] {
    let train_end = target_count(indices.len(), ratios.train).min(indices.len());
    let validation_end =
        (train_end + target_count(indices.len(), ratios.validation)).min(indices.len());
    [
        indices[..train_end].to_vec(),
        indices[train_end..validation_end].to_vec(),
        indices[validation_end..].to_vec(),
    ]
}

fn target_count(total: usize, ratio: f64) -> usize {
    f64_to_usize((usize_to_f64(total) * ratio).round())
}

fn assign_groups(mut groups: Vec<Vec<usize>>, ratios: SplitRatios) -> [Vec<usize>; 3] {
    groups.sort_by(|left, right| right.len().cmp(&left.len()).then(left[0].cmp(&right[0])));
    let targets = ratios.values();
    let total = usize_to_f64(groups.iter().map(Vec::len).sum());
    let mut output: [Vec<usize>; 3] = std::array::from_fn(|_| Vec::new());
    for group in groups {
        let partition = (0..3)
            .filter(|index| targets[*index] > 0.0)
            .min_by(|left, right| {
                let left_fill = usize_to_f64(output[*left].len()) / (total * targets[*left]);
                let right_fill = usize_to_f64(output[*right].len()) / (total * targets[*right]);
                left_fill.total_cmp(&right_fill).then(left.cmp(right))
            })
            .map_or(0, |index| index);
        output[partition].extend(group);
    }
    for partition in &mut output {
        partition.sort_unstable();
    }
    output
}

struct DisjointSet(Vec<usize>);

impl DisjointSet {
    fn new(size: usize) -> Self {
        Self((0..size).collect())
    }

    fn find(&mut self, index: usize) -> usize {
        if self.0[index] != index {
            self.0[index] = self.find(self.0[index]);
        }
        self.0[index]
    }

    fn union(&mut self, left: usize, right: usize) {
        let left = self.find(left);
        let right = self.find(right);
        if left != right {
            let (root, child) = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            self.0[child] = root;
        }
    }
}

struct SplitMix64(u64);

/// Weyl-sequence increment specified by `SplitMix64`.
const SPLITMIX_INCREMENT: u64 = 0x9E37_79B9_7F4A_7C15;
/// First avalanche multiplier specified by `SplitMix64`.
const SPLITMIX_MIX_FIRST: u64 = 0xBF58_476D_1CE4_E5B9;
/// Second avalanche multiplier specified by `SplitMix64`.
const SPLITMIX_MIX_SECOND: u64 = 0x94D0_49BB_1331_11EB;

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(SPLITMIX_INCREMENT);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(SPLITMIX_MIX_FIRST);
        value = (value ^ (value >> 27)).wrapping_mul(SPLITMIX_MIX_SECOND);
        value ^ (value >> 31)
    }
}
