//! Coordinate measurements used by frozen external-benchmark profiles.
//!
//! Correspondence and atom selection remain explicit inputs. Each RMSD is one
//! linear fit followed by one linear measurement pass. The clash distance is a
//! complete pairwise search because benchmark ligands and catalytic motifs are
//! small; callers with larger sets should supply a spatially planned distance.

use pdbiox_geom::{SuperposeError, rmsd, superpose};
use std::collections::BTreeMap;

/// Measurements used by the first `MotifBench` protocol.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotifBenchMeasurements {
    /// C-alpha self-consistency RMSD over the full scaffold.
    pub rmsd: f64,
    /// N, C-alpha, C and O RMSD over the motif.
    pub motif_rmsd: f64,
}

impl MotifBenchMeasurements {
    /// Named metrics consumed by [`crate::motifbench_1_0`].
    #[must_use]
    pub fn metrics(self) -> BTreeMap<Box<str>, f64> {
        BTreeMap::from([
            (Box::<str>::from("rmsd"), self.rmsd),
            (Box::<str>::from("motif_rmsd"), self.motif_rmsd),
        ])
    }
}

/// Measures one predicted sequence under the `MotifBench` conventions.
///
/// The full-scaffold C-alpha coordinates and motif N/C-alpha/C/O coordinates
/// are independently fitted and measured. Coordinate arrays must already be in
/// benchmark correspondence order.
///
/// # Errors
///
/// Returns the geometric reason either correspondence cannot define a unique
/// rigid fit.
pub fn measure_motifbench(
    generated_c_alpha: &[[f32; 3]],
    predicted_c_alpha: &[[f32; 3]],
    reference_motif_backbone: &[[f32; 3]],
    predicted_motif_backbone: &[[f32; 3]],
) -> Result<MotifBenchMeasurements, SuperposeError> {
    Ok(MotifBenchMeasurements {
        rmsd: fitted_rmsd(generated_c_alpha, predicted_c_alpha)?,
        motif_rmsd: fitted_rmsd(predicted_motif_backbone, reference_motif_backbone)?,
    })
}

/// Measurements used by the Atomic Motif Enzyme benchmark.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AmeMeasurements {
    /// Catalytic-residue heavy-atom RMSD after catalytic N/C-alpha/C fitting.
    pub catalytic_heavy_atom_rmsd: f64,
    /// Smallest ligand-to-predicted-backbone atom distance.
    pub ligand_backbone_min_distance: f64,
}

impl AmeMeasurements {
    /// Named metrics consumed by [`crate::ame_heavy_atom_1_0`].
    #[must_use]
    pub fn metrics(self) -> BTreeMap<Box<str>, f64> {
        BTreeMap::from([
            (
                Box::<str>::from("catalytic_heavy_atom_rmsd"),
                self.catalytic_heavy_atom_rmsd,
            ),
            (
                Box::<str>::from("ligand_backbone_min_distance"),
                self.ligand_backbone_min_distance,
            ),
        ])
    }
}

/// Why AME coordinates could not produce a benchmark measurement.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum AmeMeasurementError {
    /// Catalytic correspondence cannot define a unique rigid fit.
    #[error("catalytic alignment failed: {0:?}")]
    Alignment(SuperposeError),
    /// Heavy-atom correspondence differs in length.
    #[error("catalytic heavy-atom correspondence has different lengths")]
    HeavyAtomLengthMismatch,
    /// No ligand atom was supplied.
    #[error("AME clash measurement requires at least one ligand atom")]
    EmptyLigand,
    /// No predicted backbone atom was supplied.
    #[error("AME clash measurement requires at least one backbone atom")]
    EmptyBackbone,
}

/// Measures one predicted candidate under the AME conventions.
///
/// The transform is obtained only from catalytic N/C-alpha/C atoms, then
/// applied to the generated catalytic heavy atoms before RMSD measurement. The
/// clash observable is the minimum distance between ligand and predicted
/// backbone atoms.
///
/// # Errors
///
/// Returns an explicit error for an invalid correspondence, degenerate fit, or
/// an empty atom set that would make the clash criterion undefined.
pub fn measure_ame(
    generated_catalytic_backbone: &[[f32; 3]],
    predicted_catalytic_backbone: &[[f32; 3]],
    generated_catalytic_heavy: &[[f32; 3]],
    predicted_catalytic_heavy: &[[f32; 3]],
    predicted_ligand: &[[f32; 3]],
    predicted_backbone: &[[f32; 3]],
) -> Result<AmeMeasurements, AmeMeasurementError> {
    if generated_catalytic_heavy.len() != predicted_catalytic_heavy.len() {
        return Err(AmeMeasurementError::HeavyAtomLengthMismatch);
    }
    if predicted_ligand.is_empty() {
        return Err(AmeMeasurementError::EmptyLigand);
    }
    if predicted_backbone.is_empty() {
        return Err(AmeMeasurementError::EmptyBackbone);
    }
    let fit = superpose(generated_catalytic_backbone, predicted_catalytic_backbone)
        .map_err(AmeMeasurementError::Alignment)?;
    let mut aligned_heavy = generated_catalytic_heavy.to_vec();
    fit.transform.apply_all(&mut aligned_heavy);
    let catalytic_heavy_atom_rmsd = rmsd(&aligned_heavy, predicted_catalytic_heavy)
        .map_err(|_| AmeMeasurementError::HeavyAtomLengthMismatch)?;
    Ok(AmeMeasurements {
        catalytic_heavy_atom_rmsd,
        ligand_backbone_min_distance: minimum_distance(predicted_ligand, predicted_backbone),
    })
}

fn fitted_rmsd(mobile: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, SuperposeError> {
    superpose(mobile, reference).map(|fit| fit.rmsd)
}

fn minimum_distance(left: &[[f32; 3]], right: &[[f32; 3]]) -> f64 {
    let mut minimum_squared = f64::INFINITY;
    for &left_position in left {
        for &right_position in right {
            minimum_squared =
                minimum_squared.min(pdbiox_geom::distance_squared(left_position, right_position));
        }
    }
    minimum_squared.sqrt()
}

#[cfg(test)]
#[path = "benchmark_tests.rs"]
mod tests;
