//! Registration for the public functional-geometry facade.

use crate::fx::{
    PyAlignedMotif, PyAlignmentKind, PyAmeMeasurements, PyAtomSite, PyCandidateAggregation,
    PyComparison, PyCompatibilityProfile, PyCompatibilityVerdict, PyComponentRole, PyComponentSpec,
    PyConstraint, PyConstraintMeasurement, PyEvaluation, PyEvaluationReport,
    PyEvaluationSpecification, PyIndeterminateReason, PyMappedMotif, PyMappingSet,
    PyMeasurementOptions, PyMeasurementSet, PyMeasurementValue, PyMissingVerdict, PyMotif,
    PyMotifBenchMeasurements, PyNamedConstraint, PyProfileAlignment, PyProfileAtomSet,
    PyProfileMetric, PyRuleOutcome, PyVerdict, PyVerdictProfile, PyVerdictRule, PyVerdictStatus,
    align_intrinsic, align_with_transform, ame_heavy_atom_1_0, evaluate_motif, map_motif,
    measure_ame, measure_constraints, measure_motifbench, motifbench_1_0,
    read_evaluation_specification,
};
use pyo3::prelude::*;

pub(super) fn register_fx(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register_errors(module)?;
    module.add_class::<PyAtomSite>()?;
    module.add_class::<PyComponentRole>()?;
    module.add_class::<PyComponentSpec>()?;
    module.add_class::<PyConstraint>()?;
    module.add_class::<PyNamedConstraint>()?;
    module.add_class::<PyMotif>()?;
    module.add_class::<PyAlignmentKind>()?;
    module.add_class::<PyAlignedMotif>()?;
    module.add_class::<PyComparison>()?;
    module.add_class::<PyMissingVerdict>()?;
    module.add_class::<PyVerdictRule>()?;
    module.add_class::<PyVerdictProfile>()?;
    module.add_class::<PyRuleOutcome>()?;
    module.add_class::<PyVerdictStatus>()?;
    module.add_class::<PyVerdict>()?;
    module.add_class::<PyMappedMotif>()?;
    module.add_class::<PyMappingSet>()?;
    module.add_class::<PyMeasurementValue>()?;
    module.add_class::<PyIndeterminateReason>()?;
    module.add_class::<PyConstraintMeasurement>()?;
    module.add_class::<PyMeasurementSet>()?;
    module.add_class::<PyMeasurementOptions>()?;
    module.add_class::<PyEvaluation>()?;
    module.add_class::<PyEvaluationReport>()?;
    module.add_class::<PyProfileAtomSet>()?;
    module.add_class::<PyProfileAlignment>()?;
    module.add_class::<PyProfileMetric>()?;
    module.add_class::<PyCandidateAggregation>()?;
    module.add_class::<PyCompatibilityProfile>()?;
    module.add_class::<PyCompatibilityVerdict>()?;
    module.add_class::<PyMotifBenchMeasurements>()?;
    module.add_class::<PyAmeMeasurements>()?;
    module.add_class::<PyEvaluationSpecification>()?;

    module.add_function(wrap_pyfunction!(align_intrinsic, module)?)?;
    module.add_function(wrap_pyfunction!(align_with_transform, module)?)?;
    module.add_function(wrap_pyfunction!(map_motif, module)?)?;
    module.add_function(wrap_pyfunction!(measure_constraints, module)?)?;
    module.add_function(wrap_pyfunction!(evaluate_motif, module)?)?;
    module.add_function(wrap_pyfunction!(measure_motifbench, module)?)?;
    module.add_function(wrap_pyfunction!(measure_ame, module)?)?;
    module.add_function(wrap_pyfunction!(motifbench_1_0, module)?)?;
    module.add_function(wrap_pyfunction!(ame_heavy_atom_1_0, module)?)?;
    module.add_function(wrap_pyfunction!(read_evaluation_specification, module)?)?;
    Ok(())
}

fn register_errors(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::fx::register(module)
}
