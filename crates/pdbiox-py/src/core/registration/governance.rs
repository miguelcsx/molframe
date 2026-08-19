//! Registration for ready-to-run governed analysis projections.

use crate::analysis::{
    PyAltlocOccupancyIssue, PyAltlocOccupancyOptions, PyAltlocOccupancyRecord,
    PyAltlocOccupancyReport, PyBFactorDistribution, PyBFactorOutlier, PyBaseFrame,
    PyCcdCompletenessReport, PyDielectricOptions, PyDielectricResult, PyHelicalOptions,
    PyHelicalParameters, PyPlaneRestraint, PyPlaneRestraintFlag, PyPlaneRestraintReport,
    PyResidueAtomCompleteness, PyTlsBFactorFlag, PyTlsBFactorReport, PyTlsGroup, PyTlsModel,
    PyWaterDynamicsOptions, PyWaterLag, analyse_altloc_occupancy, analyse_b_factor_distribution,
    analyse_ccd_completeness, analyse_centre_of_mass_radial_distribution,
    analyse_dielectric_from_dipoles, analyse_helical_parameters, analyse_helical_steps,
    analyse_plane_restraints, analyse_tls_b_factor_consistency, analyse_water_dynamics,
    b_factor_distribution, dielectric_from_dipoles, helical_parameters, helical_steps,
    tls_b_factor_consistency, validate_altloc_occupancy, validate_ccd_completeness,
    validate_plane_restraints, water_dynamics,
};
use pyo3::prelude::*;

pub(super) fn register_governance(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyWaterDynamicsOptions>()?;
    module.add_class::<PyWaterLag>()?;
    module.add_class::<PyDielectricOptions>()?;
    module.add_class::<PyDielectricResult>()?;
    module.add_class::<PyBaseFrame>()?;
    module.add_class::<PyHelicalOptions>()?;
    module.add_class::<PyHelicalParameters>()?;
    module.add_class::<PyAltlocOccupancyOptions>()?;
    module.add_class::<PyAltlocOccupancyIssue>()?;
    module.add_class::<PyAltlocOccupancyRecord>()?;
    module.add_class::<PyAltlocOccupancyReport>()?;
    module.add_class::<PyResidueAtomCompleteness>()?;
    module.add_class::<PyCcdCompletenessReport>()?;
    module.add_class::<PyPlaneRestraint>()?;
    module.add_class::<PyPlaneRestraintFlag>()?;
    module.add_class::<PyPlaneRestraintReport>()?;
    module.add_class::<PyBFactorOutlier>()?;
    module.add_class::<PyBFactorDistribution>()?;
    module.add_class::<PyTlsModel>()?;
    module.add_class::<PyTlsGroup>()?;
    module.add_class::<PyTlsBFactorFlag>()?;
    module.add_class::<PyTlsBFactorReport>()?;
    module.add_function(wrap_pyfunction!(water_dynamics, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_water_dynamics, module)?)?;
    module.add_function(wrap_pyfunction!(dielectric_from_dipoles, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_dielectric_from_dipoles, module)?)?;
    module.add_function(wrap_pyfunction!(helical_parameters, module)?)?;
    module.add_function(wrap_pyfunction!(helical_steps, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_helical_parameters, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_helical_steps, module)?)?;
    module.add_function(wrap_pyfunction!(
        analyse_centre_of_mass_radial_distribution,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(validate_altloc_occupancy, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_altloc_occupancy, module)?)?;
    module.add_function(wrap_pyfunction!(validate_ccd_completeness, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_ccd_completeness, module)?)?;
    module.add_function(wrap_pyfunction!(validate_plane_restraints, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_plane_restraints, module)?)?;
    module.add_function(wrap_pyfunction!(b_factor_distribution, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_b_factor_distribution, module)?)?;
    module.add_function(wrap_pyfunction!(tls_b_factor_consistency, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_tls_b_factor_consistency, module)?)?;
    Ok(())
}
