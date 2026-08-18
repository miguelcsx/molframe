/// Agglomerative inter-cluster distance definition.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Linkage {
    /// Minimum pair distance.
    Single,
    /// Maximum pair distance.
    Complete,
    /// Arithmetic mean pair distance.
    Average,
}

/// Deterministic cluster membership and medoid representatives.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Clustering {
    /// Cluster label per observation; `None` denotes DBSCAN noise.
    pub labels: Vec<Option<usize>>,
    /// Sorted observation indices for each cluster.
    pub members: Vec<Vec<usize>>,
    /// Minimum-total-distance member of each cluster.
    pub medoids: Vec<usize>,
}

/// Packed mutable distances between active agglomerative clusters.
///
/// Only the upper triangle is stored because the validated source matrix is
/// symmetric and cluster-linkage updates preserve symmetry.
struct PairDistances {
    size: usize,
    offsets: Vec<usize>,
    values: Vec<f64>,
}

impl PairDistances {
    /// Copies the source matrix's upper triangle into mutable packed storage.
    fn from_matrix(distances: &EnsembleDistanceMatrix) -> Result<Self, EnsembleGeometryError> {
        let size = distances.size;
        let pair_count = checked_pair_count(size).ok_or(EnsembleGeometryError::InvalidParameter)?;

        let mut offsets = Vec::new();
        offsets
            .try_reserve_exact(size)
            .map_err(|_| EnsembleGeometryError::InvalidParameter)?;

        let mut values = Vec::new();
        values
            .try_reserve_exact(pair_count)
            .map_err(|_| EnsembleGeometryError::InvalidParameter)?;

        for left in 0..size {
            offsets.push(values.len());
            let row = matrix_row(distances, left)?;

            for right in left + 1..size {
                let value = row
                    .get(right)
                    .copied()
                    .ok_or(EnsembleGeometryError::DimensionMismatch)?;
                values.push(canonical_distance(value));
            }
        }

        Ok(Self {
            size,
            offsets,
            values,
        })
    }

    /// Returns the current distance between two distinct cluster identifiers.
    fn get(&self, left: usize, right: usize) -> Result<f64, EnsembleGeometryError> {
        let index = self
            .index(left, right)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        self.values
            .get(index)
            .copied()
            .ok_or(EnsembleGeometryError::InvalidParameter)
    }

    /// Replaces the current distance between two distinct cluster identifiers.
    fn set(&mut self, left: usize, right: usize, value: f64) -> Result<(), EnsembleGeometryError> {
        let index = self
            .index(left, right)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        let slot = self
            .values
            .get_mut(index)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        *slot = canonical_distance(value);
        Ok(())
    }

    /// Resolves an unordered cluster pair into packed upper-triangle storage.
    fn index(&self, left: usize, right: usize) -> Option<usize> {
        if left == right || left >= self.size || right >= self.size {
            return None;
        }

        let (left, right) = ordered_pair(left, right);
        let offset = *self.offsets.get(left)?;
        let delta = right.checked_sub(left)?.checked_sub(1)?;

        offset.checked_add(delta)
    }
}

/// One heap entry representing a candidate agglomerative merge.
///
/// Version fields make candidates involving an already-merged cluster lazily
/// detectable without deleting arbitrary entries from the binary heap.
#[derive(Clone, Copy, Debug)]
struct MergeCandidate {
    distance: f64,
    left: usize,
    right: usize,
    left_version: usize,
    right_version: usize,
}

impl MergeCandidate {
    /// Creates a candidate with a canonical zero representation.
    fn new(
        distance: f64,
        left: usize,
        right: usize,
        left_version: usize,
        right_version: usize,
    ) -> Self {
        Self {
            distance: canonical_distance(distance),
            left,
            right,
            left_version,
            right_version,
        }
    }
}

impl PartialEq for MergeCandidate {
    /// Compares the complete heap identity of two merge candidates.
    fn eq(&self, other: &Self) -> bool {
        self.distance.to_bits() == other.distance.to_bits()
            && self.left == other.left
            && self.right == other.right
            && self.left_version == other.left_version
            && self.right_version == other.right_version
    }
}

impl Eq for MergeCandidate {}

impl PartialOrd for MergeCandidate {
    /// Delegates candidate ordering to the total heap ordering.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MergeCandidate {
    /// Orders the min-distance, lexicographically earliest pair first.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.left.cmp(&self.left))
            .then_with(|| other.right.cmp(&self.right))
            .then_with(|| other.left_version.cmp(&self.left_version))
            .then_with(|| other.right_version.cmp(&self.right_version))
    }
}
