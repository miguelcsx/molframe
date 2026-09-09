//! Crystallographic Python projections split by ownership and operation type.

#[path = "assembly.rs"]
mod assembly;
#[path = "functions.rs"]
mod functions;
#[path = "ncs.rs"]
mod ncs;
#[path = "reflection.rs"]
mod reflection;
#[path = "reflection_metadata.rs"]
mod reflection_metadata;
#[path = "types.rs"]
mod types;

pub(crate) use reflection::PyReflectionTable;

use pyo3::prelude::*;

pub(crate) use assembly::{
    PyAssemblyNeighbor, PyAssemblySet, PyAssemblyView, PyAtomInstance, PyChainInstance,
};
pub(crate) use functions::{collect_crystal_neighbors, lower_assemblies, lower_symmetry};
pub(crate) use ncs::{PyCrystalNeighbor, PyNcsCode, PyNcsOperator, PyNcsSet, PyNcsView, lower_ncs};
pub(crate) use types::{
    PyAffineTransform, PyAssemblyDef, PyGenerator, PyOperExpression, PyOperator, PyRational,
    PySymmetrySet,
};

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "DEFAULT_CRYSTAL_IMAGE_LIMIT",
        pdbiox::DEFAULT_CRYSTAL_IMAGE_LIMIT,
    )?;
    module.add("DEFAULT_INSTANCE_LIMIT", pdbiox::DEFAULT_INSTANCE_LIMIT)?;
    module.add_class::<PyAffineTransform>()?;
    module.add_class::<PyRational>()?;
    module.add_class::<PyOperExpression>()?;
    module.add_class::<PySymmetrySet>()?;
    module.add_class::<PyOperator>()?;
    module.add_class::<PyGenerator>()?;
    module.add_class::<PyAssemblyDef>()?;
    module.add_class::<PyAssemblySet>()?;
    module.add_class::<PyChainInstance>()?;
    module.add_class::<PyAtomInstance>()?;
    module.add_class::<PyAssemblyNeighbor>()?;
    module.add_class::<PyAssemblyView>()?;
    module.add_class::<PyNcsCode>()?;
    module.add_class::<PyNcsOperator>()?;
    module.add_class::<PyNcsSet>()?;
    module.add_class::<PyNcsView>()?;
    module.add_class::<PyCrystalNeighbor>()?;
    module.add_function(wrap_pyfunction!(lower_assemblies, module)?)?;
    module.add_function(wrap_pyfunction!(lower_ncs, module)?)?;
    module.add_function(wrap_pyfunction!(lower_symmetry, module)?)?;
    module.add_function(wrap_pyfunction!(collect_crystal_neighbors, module)?)?;
    reflection::register(module)?;
    super::grids::register(module)?;
    super::bricks::register(module)?;
    Ok(())
}
