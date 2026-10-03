//! Explicit bond inference for Python structures.

use super::PyStructure;
use pyo3::prelude::*;

/// Infers bonds by covalent-radius distance, keeping every existing bond.
///
/// Files such as mmCIF usually list only the connections between residues
/// and ligands, so a structure read from one draws no bonds inside a residue
/// until they are inferred.
pub(super) fn infer(
    py: Python<'_>,
    structure: &PyStructure,
    scale: f32,
    lower_bound: f32,
    across_chains: bool,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyStructure> {
    let options = molframe::BondInference {
        scale,
        lower_bound,
        exclude_across_chains: !across_chains,
        ..molframe::BondInference::default()
    };
    let source = structure.inner.clone();
    crate::execution::run(py, context, move |context| {
        molframe::infer_bonds(&source, options, context)
    })?
    .map(|report| PyStructure {
        inner: report.structure,
    })
    .map_err(|findings| super::findings_error(&findings))
}
