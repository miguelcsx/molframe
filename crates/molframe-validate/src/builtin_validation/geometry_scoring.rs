//! Bond-length and bond-angle z-scores against the reference targets.
//!
//! Every standard residue is measured against its own targets and against the
//! peptide-link terms to the previous residue, each as `(observed - ideal) /
//! sigma`. The per-residue arrays keep the worst signed z of each kind, the
//! outlier list names the atoms behind every z beyond the cutoff, and the
//! root-mean-square z over all terms summarises the model.
//!
//! The targets come from the resolved geometry table, whose identity the result
//! records.
//!
//! Complexity: O(residues + terms), with each term a constant-time atom lookup
//! and one distance or angle. Parallel blocks of 256 residues; no allocation
//! per residue, since the peptide-link targets are iterated lazily.

use molframe_chem::TableIdentity;
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};

use super::reference_tables::{GeometryTargets, Target};
use super::{BuiltinError, Layout, ReferenceSet};
use crate::numeric::{f64_to_f32, usize_to_f64};

/// Residues per parallel block.
const BLOCK: usize = 256;
/// Atom slot value for the unused third atom of a bond.
pub const NO_ATOM: u32 = u32::MAX;

/// One bond or angle beyond the z cutoff.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeometryOutlier {
    /// Residue the term belongs to; for a link term, the later residue.
    pub residue: u32,
    /// Atom indices in bonded order; the third is [`NO_ATOM`] for a bond.
    pub atoms: [u32; 3],
    /// Observed length (ångström) or angle (degrees).
    pub observed: f32,
    /// Ideal value.
    pub expected: f32,
    /// Signed z-score.
    pub z: f32,
}

impl GeometryOutlier {
    /// Whether this is an angle rather than a bond length.
    #[must_use]
    pub const fn is_angle(&self) -> bool {
        self.atoms[2] != NO_ATOM
    }
}

/// Bond and angle z-scores for a whole structure.
#[derive(Clone, Debug, Default)]
pub struct GeometryScores {
    /// Signed bond z with the largest magnitude per residue, `NaN` if none.
    pub worst_bond_z: Vec<f32>,
    /// Signed angle z with the largest magnitude per residue, `NaN` if none.
    pub worst_angle_z: Vec<f32>,
    /// Terms beyond the cutoff, ordered by residue.
    pub outliers: Vec<GeometryOutlier>,
    /// Number of bonds measured.
    pub bond_count: usize,
    /// Number of angles measured.
    pub angle_count: usize,
    /// Root-mean-square z of the bonds, `NaN` if none.
    pub bond_rmsz: f64,
    /// Root-mean-square z of the angles, `NaN` if none.
    pub angle_rmsz: f64,
    /// The target table the z-scores were computed against.
    pub tables: Vec<TableIdentity>,
}

#[derive(Default)]
struct Block {
    bond: Vec<f32>,
    angle: Vec<f32>,
    outliers: Vec<GeometryOutlier>,
    counts: [usize; 2],
    squares: [f64; 2],
}

fn locate(
    layout: &Layout<'_>,
    target: &Target,
    residue: usize,
    link: bool,
) -> Option<[(u32, [f32; 3]); 3]> {
    let mut found = [(NO_ATOM, [0.0f32; 3]); 3];
    let mut seen_nitrogen = false;
    for (slot, name) in found.iter_mut().zip(&target.names) {
        seen_nitrogen |= name == "N";
        let owner = if link && !seen_nitrogen {
            residue.checked_sub(1)?
        } else {
            residue
        };
        *slot = layout.atom(owner, name)?;
    }
    Some(found)
}

fn measure(target: &Target, found: &[(u32, [f32; 3]); 3]) -> Option<f64> {
    let [a, b, c] = found;
    if target.is_angle() {
        molframe_geom::angle(a.1, b.1, c.1).map(f64::to_degrees)
    } else {
        Some(molframe_geom::distance(a.1, b.1))
    }
}

fn score_residue(
    layout: &Layout<'_>,
    targets: &GeometryTargets,
    residue: usize,
    cutoff: f64,
    block: &mut Block,
) {
    let mut worst = [0.0f64; 2];
    let name = layout.name(residue);
    let linked = layout.peptide_linked(residue);
    let own = targets.residue(name).iter().map(|target| (target, false));
    let link = targets
        .link(name)
        .filter(|_| linked)
        .map(|target| (target, true));
    for (target, is_link) in own.chain(link) {
        let Some(found) = locate(layout, target, residue, is_link) else {
            continue;
        };
        let Some(observed) = measure(target, &found) else {
            continue;
        };
        let z = (observed - target.mean) / target.sigma;
        let kind = usize::from(target.is_angle());
        block.counts[kind] += 1;
        block.squares[kind] += z * z;
        if z.abs() > worst[kind].abs() {
            worst[kind] = z;
        }
        if z.abs() > cutoff {
            block.outliers.push(GeometryOutlier {
                residue: match u32::try_from(residue) {
                    Ok(index) => index,
                    Err(_) => u32::MAX,
                },
                atoms: [
                    found[0].0,
                    found[1].0,
                    if target.is_angle() {
                        found[2].0
                    } else {
                        NO_ATOM
                    },
                ],
                observed: f64_to_f32(observed),
                expected: f64_to_f32(target.mean),
                z: f64_to_f32(z),
            });
        }
    }
    let present = |kind: usize| {
        if block.counts[kind] > 0 && worst[kind] != 0.0 {
            f64_to_f32(worst[kind])
        } else {
            f32::NAN
        }
    };
    block.bond.push(present(0));
    block.angle.push(present(1));
}

/// Scores every bond and angle that has a target.
///
/// Targets come from `reference`, which the result records. `cutoff` is the
/// z magnitude above which a term is listed as an outlier; 4 is the
/// conventional choice. Runs in `O(residues)` time across the
/// context's workers; block results are joined in order so the output is
/// identical at any worker count.
///
/// # Errors
///
/// Returns [`BuiltinError`] if a worker panics.
pub fn geometry_scores(
    layout: &Layout<'_>,
    reference: &ReferenceSet,
    cutoff: f64,
    context: &ExecutionContext,
) -> Result<GeometryScores, BuiltinError> {
    let targets = reference.geometry()?;
    let plan = BlockPlan::new(layout.len(), BLOCK);
    let blocks = map_blocks_in(plan, context, |_, range| {
        let mut block = Block::default();
        for residue in range {
            // The worst-z bookkeeping above is per residue, so a residue
            // without terms still pushes its NaN pair.
            score_residue(layout, &targets.table, residue, cutoff, &mut block);
        }
        block
    })?;
    let mut scores = GeometryScores {
        tables: vec![targets.identity.clone()],
        ..GeometryScores::default()
    };
    let mut squares = [0.0f64; 2];
    for block in blocks {
        scores.worst_bond_z.extend(block.bond);
        scores.worst_angle_z.extend(block.angle);
        scores.outliers.extend(block.outliers);
        scores.bond_count += block.counts[0];
        scores.angle_count += block.counts[1];
        squares[0] += block.squares[0];
        squares[1] += block.squares[1];
    }
    let rms = |sum: f64, count: usize| {
        if count == 0 {
            f64::NAN
        } else {
            (sum / usize_to_f64(count)).sqrt()
        }
    };
    scores.bond_rmsz = rms(squares[0], scores.bond_count);
    scores.angle_rmsz = rms(squares[1], scores.angle_count);
    Ok(scores)
}

#[cfg(test)]
#[path = "geometry_scoring_tests.rs"]
mod tests;
