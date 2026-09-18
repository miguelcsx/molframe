//! Registration for internal-coordinate projections.

use crate::internal_coordinates::{
    PyBatFrame, PyDihedron, PyHedron, PyInternalAtom, PyInternalCoordinates, internal_coordinates,
    place_atom,
};
use pyo3::prelude::*;

pub(super) fn register_ic(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyHedron>()?;
    module.add_class::<PyDihedron>()?;
    module.add_class::<PyInternalAtom>()?;
    module.add_class::<PyBatFrame>()?;
    module.add_class::<PyInternalCoordinates>()?;
    module.add_function(wrap_pyfunction!(internal_coordinates, module)?)?;
    module.add_function(wrap_pyfunction!(place_atom, module)?)?;
    Ok(())
}
