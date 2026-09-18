//! Registration for geometry result classes.

use crate::geometry::{
    PyAxes, PyBackboneCoordinates, PyBackboneFrame, PyBackboneTorsions, PyCircularSummary,
    PyEigenOptions, PyHelixGeometry, PyPlane,
};
use pyo3::prelude::*;

pub(super) fn register_geometry_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAxes>()?;
    module.add_class::<PyEigenOptions>()?;
    module.add_class::<PyPlane>()?;
    module.add_class::<PyCircularSummary>()?;
    module.add_class::<PyBackboneCoordinates>()?;
    module.add_class::<PyBackboneFrame>()?;
    module.add_class::<PyBackboneTorsions>()?;
    module.add_class::<PyHelixGeometry>()?;
    crate::geometry::register(module)?;
    Ok(())
}
