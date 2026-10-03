//! Mechanical adapters for reciprocal geometry and exact reflection symmetry.

use crate::bindings::PyStructure;
use molframe::UnitCell;
use molframe::crystal::{AssemblyExt as _, AssemblyView};
use molframe::crystal::{CellTransform, ReflectionSymmetry, SymmetrySet, space_group_setting};
use numpy::{Complex64, PyArray1, PyArray2, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
use pyo3::prelude::*;

#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "UnitCell",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PyUnitCell(pub(crate) CellTransform);

#[pymethods]
impl PyUnitCell {
    #[new]
    fn new(lengths: [f64; 3], angles: [f64; 3]) -> PyResult<Self> {
        CellTransform::new(&UnitCell { lengths, angles })
            .map(Self)
            .map_err(crate::error::kernel)
    }

    fn reciprocal_vector(&self, hkl: [i32; 3]) -> [f64; 3] {
        self.0.reciprocal_vector(hkl)
    }

    fn reciprocal_spacing_squared(&self, hkl: [i32; 3]) -> f64 {
        self.0.reciprocal_spacing_squared(hkl)
    }

    fn d_spacing(&self, hkl: [i32; 3]) -> f64 {
        self.0.d_spacing(hkl)
    }
}

#[derive(Clone, Debug)]
#[pyclass(
    name = "SpaceGroup",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PySpaceGroup(pub(crate) SymmetrySet);

#[pymethods]
impl PySpaceGroup {
    #[new]
    fn new(hall_number: u16) -> PyResult<Self> {
        space_group_setting(hall_number)
            .map(|setting| Self(setting.symmetry_set()))
            .map_err(crate::error::kernel)
    }

    #[getter]
    fn hall_symbol(&self) -> Option<&str> {
        self.0.hall.as_deref()
    }

    #[getter]
    fn hall_number(&self) -> Option<u16> {
        self.0.hall_number
    }

    fn reflection_symmetry(&self, hkl: [i32; 3]) -> PyResult<PyReflectionSymmetry> {
        self.0
            .reflection_symmetry(hkl)
            .map(PyReflectionSymmetry)
            .map_err(crate::error::kernel)
    }
}

#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "ReflectionSymmetry",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PyReflectionSymmetry(ReflectionSymmetry);

#[pymethods]
impl PyReflectionSymmetry {
    #[getter]
    fn centric(&self) -> bool {
        self.0.centric
    }
    #[getter]
    fn systematically_absent(&self) -> bool {
        self.0.systematically_absent
    }
    #[getter]
    fn epsilon_factor(&self) -> usize {
        self.0.epsilon_factor
    }
}

/// X-ray structure factors of the first model, one complex value per `(h, k, l)` row.
///
/// The cell, space group, occupancies and displacement parameters come from the
/// structure itself; the model is taken to be the asymmetric unit.
#[pyfunction]
fn structure_factors<'py>(
    py: Python<'py>,
    structure: &PyStructure,
    hkl: &Bound<'py, PyArray2<i32>>,
) -> PyResult<Bound<'py, PyArray1<Complex64>>> {
    let hkl = hkl.readonly();
    let shape = hkl.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(crate::error::value("hkl must have shape (n, 3)"));
    }
    let rows: Vec<[i32; 3]> = hkl
        .as_array()
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    let structure = structure.inner.clone();
    let values = py
        .detach(move || molframe::crystal::structure_factors(structure.engine(), &rows))
        .map_err(crate::error::kernel)?;
    let values: Vec<Complex64> = values
        .into_iter()
        .map(|value| Complex64::new(value.re, value.im))
        .collect();
    Ok(values.to_pyarray(py))
}

/// One placement of an assembly: a transform and the chains it applies to.
#[derive(Clone, Debug)]
#[pyclass(
    name = "AssemblyInstance",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PyAssemblyInstance {
    matrix: [f64; 16],
    chains: Vec<String>,
}

#[pymethods]
impl PyAssemblyInstance {
    /// Column-major 4×4 affine matrix in ångström.
    #[getter]
    fn matrix(&self) -> Vec<f64> {
        self.matrix.to_vec()
    }

    /// `label_asym_id` of every chain the transform applies to.
    #[getter]
    fn chains(&self) -> Vec<String> {
        self.chains.clone()
    }

    fn __repr__(&self) -> String {
        format!("AssemblyInstance(chains={})", self.chains.join(","))
    }
}

/// Identifiers of the biological assemblies a structure declares.
#[pyfunction]
fn assemblies(structure: &PyStructure) -> Vec<String> {
    match structure.inner.engine().assembly_set() {
        Some(set) => set.assemblies().map(|each| each.id.to_string()).collect(),
        // A structure that declares no assemblies has none to list.
        None => Vec::new(),
    }
}

/// The placements that make up one biological assembly.
///
/// Each carries a transform and the chains it applies to; the transform of
/// every chain is the product of the operators the entry lists for it.
#[pyfunction]
fn assembly(structure: &PyStructure, id: &str) -> PyResult<Vec<PyAssemblyInstance>> {
    let engine = structure.inner.engine();
    let set = engine
        .assembly_set()
        .ok_or_else(|| crate::error::value("the structure declares no biological assemblies"))?;
    let view = AssemblyView::new(engine, set, id).map_err(crate::error::kernel)?;
    view.groups_by_transform()
        .into_iter()
        .map(|(transform, chains)| {
            let mut matrix = [0.0_f64; 16];
            for column in 0..3 {
                for row in 0..3 {
                    matrix[column * 4 + row] = transform.rotation[row][column];
                }
            }
            matrix[12..15].copy_from_slice(&transform.translation);
            matrix[15] = 1.0;
            let labels = chains
                .into_iter()
                .map(|chain| {
                    engine
                        .chain(chain)
                        .and_then(molframe::ChainRef::label)
                        .map(str::to_owned)
                        .ok_or_else(|| crate::error::value("an assembly chain has no label"))
                })
                .collect::<PyResult<Vec<_>>>()?;
            Ok(PyAssemblyInstance {
                matrix,
                chains: labels,
            })
        })
        .collect()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::crystal_links::register(module)?;
    module.add_class::<PyAssemblyInstance>()?;
    module.add_function(wrap_pyfunction!(assemblies, module)?)?;
    module.add_function(wrap_pyfunction!(assembly, module)?)?;
    module.add_function(wrap_pyfunction!(structure_factors, module)?)?;
    module.add_class::<PyUnitCell>()?;
    module.add_class::<PySpaceGroup>()?;
    module.add_class::<PyReflectionSymmetry>()?;
    crate::crystal_statistics::register(module)?;
    crate::crystal_reduction::register(module)
}
