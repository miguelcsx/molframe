//! Registration for direct public facade projections.

use crate::cif_document::{
    PyCifCategory, PyCifColumn, PyCifDataBlock, PyCifDocument, PyCifQuoting, PyCifValue,
    read_document, write_preserving,
};
use crate::config::{
    PyApplicationConfiguration, PyChemistryConfiguration, PyOutputConfiguration, PyPolicyOverrides,
    read_configuration, read_policy, read_policy_overrides,
};
use crate::facade::{
    PyBackboneTorsionRecord, PyBondInference, PyBondInferenceReport, PyProteinAlphaTrace,
    PySideChainTorsionRecord, PySideChainTorsionReport, default_limits, infer_bonds,
    structure_backbone_torsions, structure_backbone_torsions_model, structure_protein_alpha_traces,
    structure_side_chain_torsions,
};
use crate::pdb_headers::register as register_pdb_headers;
use pyo3::prelude::*;

pub(super) fn register_facade(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::difference::register(module)?;
    crate::cif_parser::register(module)?;
    crate::cif_lower::register(module)?;
    crate::cif_pdbml::register(module)?;
    crate::cif_small::register(module)?;
    crate::cif_write::register(module)?;
    register_pdb_headers(module)?;
    module.add_class::<PyBondInference>()?;
    module.add_class::<PyBondInferenceReport>()?;
    module.add_class::<PyProteinAlphaTrace>()?;
    module.add_class::<PyBackboneTorsionRecord>()?;
    module.add_class::<PySideChainTorsionRecord>()?;
    module.add_class::<PySideChainTorsionReport>()?;
    module.add_class::<PyCifValue>()?;
    module.add_class::<PyCifQuoting>()?;
    crate::cif_rows::register(module)?;
    crate::cif_lexer::register(module)?;
    module.add_class::<PyCifColumn>()?;
    module.add_class::<PyCifCategory>()?;
    module.add_class::<PyCifDataBlock>()?;
    module.add_class::<PyCifDocument>()?;
    module.add_class::<PyPolicyOverrides>()?;
    module.add_class::<PyOutputConfiguration>()?;
    module.add_class::<PyChemistryConfiguration>()?;
    module.add_class::<PyApplicationConfiguration>()?;
    module.add_function(wrap_pyfunction!(default_limits, module)?)?;
    module.add_function(wrap_pyfunction!(infer_bonds, module)?)?;
    module.add_function(wrap_pyfunction!(structure_protein_alpha_traces, module)?)?;
    module.add_function(wrap_pyfunction!(structure_backbone_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(structure_backbone_torsions_model, module)?)?;
    module.add_function(wrap_pyfunction!(structure_side_chain_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(read_document, module)?)?;
    module.add_function(wrap_pyfunction!(write_preserving, module)?)?;
    module.add_function(wrap_pyfunction!(read_configuration, module)?)?;
    module.add_function(wrap_pyfunction!(read_policy, module)?)?;
    module.add_function(wrap_pyfunction!(read_policy_overrides, module)?)?;
    Ok(())
}
