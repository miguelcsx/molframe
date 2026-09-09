//! Mechanical projections for anisotropic networks and NMD interchange.

use crate::geometry::coordinates;
use crate::graph::PySpatialBackend;
use crate::query::PySelection;
use crate::structure::PyStructure;
use numpy::ndarray::{Array1, Array2, Array3};
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyArray3, PyReadonlyArray1, PyReadonlyArray2,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::periodic::periodic_box;

#[pyclass(name = "AnmOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAnmOptions(pub(crate) pdbiox::analysis::AnmOptions);

#[pymethods]
impl PyAnmOptions {
    #[new]
    fn new(
        contact_distance: f32,
        mode_count: usize,
        zero_mode_tolerance: f64,
        memory_limit_bytes: usize,
        backend: PySpatialBackend,
    ) -> Self {
        Self(pdbiox::analysis::AnmOptions {
            contact_distance,
            mode_count,
            zero_mode_tolerance,
            memory_limit_bytes,
            backend: backend.into(),
        })
    }
}

#[pyclass(name = "AnisotropicNetworkModel", frozen, skip_from_py_object)]
pub(crate) struct PyAnisotropicNetworkModel {
    sites: Py<PyArray1<u32>>,
    eigenvalues: Py<PyArray1<f64>>,
    modes: Py<PyArray3<f64>>,
    native: pdbiox::analysis::AnisotropicNetworkModel,
}

#[pymethods]
impl PyAnisotropicNetworkModel {
    #[getter]
    fn sites(&self, py: Python<'_>) -> Py<PyArray1<u32>> {
        self.sites.clone_ref(py)
    }
    #[getter]
    fn eigenvalues(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.eigenvalues.clone_ref(py)
    }
    #[getter]
    fn modes(&self, py: Python<'_>) -> Py<PyArray3<f64>> {
        self.modes.clone_ref(py)
    }
    #[getter]
    fn zero_modes(&self) -> usize {
        self.native.zero_modes
    }

    fn fluctuations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_vec(self.native.fluctuations()).into_pyarray(py)
    }

    fn project<'py>(
        &self,
        py: Python<'py>,
        values: PyReadonlyArray2<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(Array1::from_vec(self.native.project(&vectors_f64(values)?)).into_pyarray(py))
    }

    fn displace<'py>(
        &self,
        py: Python<'py>,
        positions: PyReadonlyArray2<'_, f32>,
        amplitudes: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let positions = coordinates(positions)?;
        let mut output = vec![[0.0; 3]; positions.len()];
        self.native
            .displace(&positions, amplitudes.as_slice()?, &mut output);
        matrix_f32(py, output)
    }
}

impl PyAnisotropicNetworkModel {
    fn new(py: Python<'_>, native: pdbiox::analysis::AnisotropicNetworkModel) -> PyResult<Self> {
        let values = native.modes.iter().flatten().flatten().copied().collect();
        let modes = Array3::from_shape_vec((native.modes.len(), native.sites.len(), 3), values)
            .map_err(value_error)?
            .into_pyarray(py)
            .unbind();
        Ok(Self {
            sites: Array1::from_vec(native.sites.clone())
                .into_pyarray(py)
                .unbind(),
            eigenvalues: Array1::from_vec(native.eigenvalues.clone())
                .into_pyarray(py)
                .unbind(),
            modes,
            native,
        })
    }
}

#[pyclass(name = "NormalMode", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNormalMode(pub(crate) pdbiox::analysis::NormalMode);

#[pymethods]
impl PyNormalMode {
    #[new]
    fn new(index: i32, scale: f64, displacements: PyReadonlyArray2<'_, f64>) -> PyResult<Self> {
        Ok(Self(pdbiox::analysis::NormalMode {
            index,
            scale,
            displacements: vectors_f64(displacements)?,
        }))
    }
    #[getter]
    fn index(&self) -> i32 {
        self.0.index
    }
    #[getter]
    fn scale(&self) -> f64 {
        self.0.scale
    }
    #[getter]
    fn displacements<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        matrix_f64(py, self.0.displacements.clone())
    }
}

#[pyclass(name = "NormalModeSet", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNormalModeSet(pub(crate) pdbiox::analysis::NormalModeSet);

#[pymethods]
impl PyNormalModeSet {
    #[getter]
    fn name(&self) -> Option<&str> {
        self.0.name.as_deref()
    }
    #[getter]
    fn atom_names(&self) -> Vec<&str> {
        self.0.atom_names.iter().map(AsRef::as_ref).collect()
    }
    #[getter]
    fn residue_names(&self) -> Vec<&str> {
        self.0.residue_names.iter().map(AsRef::as_ref).collect()
    }
    #[getter]
    fn residue_ids(&self) -> &[i32] {
        &self.0.residue_ids
    }
    #[getter]
    fn chain_ids(&self) -> Vec<&str> {
        self.0.chain_ids.iter().map(AsRef::as_ref).collect()
    }
    #[getter]
    fn b_factors(&self) -> &[f64] {
        &self.0.b_factors
    }
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        matrix_f64(py, self.0.coordinates.clone())
    }
    #[getter]
    fn modes(&self) -> Vec<PyNormalMode> {
        self.0.modes.iter().cloned().map(PyNormalMode).collect()
    }
    fn __len__(&self) -> usize {
        self.0.len()
    }
    fn displace<'py>(
        &self,
        py: Python<'py>,
        amplitudes: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let mut output = vec![[0.0; 3]; self.0.len()];
        self.0.displace(amplitudes.as_slice()?, &mut output);
        matrix_f64(py, output)
    }
}

#[pyfunction]
#[pyo3(signature = (structure, sites, options, *, periodic=false))]
pub(crate) fn anisotropic_network_model(
    py: Python<'_>,
    structure: &PyStructure,
    sites: &PySelection,
    options: PyAnmOptions,
    periodic: bool,
) -> PyResult<PyAnisotropicNetworkModel> {
    let structure = structure.structure().clone();
    let sites = sites.inner.clone();
    let model = py
        .detach(move || {
            let cell = periodic_box(&structure, periodic)?;
            pdbiox::analysis::anisotropic_network_model(
                structure.positions(),
                &sites,
                options.0,
                cell.as_ref(),
                &crate::core::execution::default_context(),
            )
            .map_err(|error| error.to_string())
        })
        .map_err(PyValueError::new_err)?;
    PyAnisotropicNetworkModel::new(py, model)
}

#[pyfunction]
pub(crate) fn read_nmd(py: Python<'_>, data: &[u8]) -> PyResult<PyNormalModeSet> {
    py.detach(move || -> PyResult<PyNormalModeSet> {
        pdbiox::analysis::read_nmd(data)
            .map(PyNormalModeSet)
            .map_err(value_error)
    })
}

#[pyfunction]
pub(crate) fn write_nmd(py: Python<'_>, set: &PyNormalModeSet) -> String {
    py.detach(move || -> String { pdbiox::analysis::write_nmd(&set.0) })
}

fn vectors_f64(array: PyReadonlyArray2<'_, f64>) -> PyResult<Vec<[f64; 3]>> {
    if array.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("vectors must have shape (n, 3)"));
    }
    Ok(array
        .as_array()
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

fn matrix_f64(py: Python<'_>, values: Vec<[f64; 3]>) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let rows = values.len();
    Array2::from_shape_vec((rows, 3), values.into_iter().flatten().collect())
        .map_err(value_error)
        .map(|value| value.into_pyarray(py))
}

fn matrix_f32(py: Python<'_>, values: Vec<[f32; 3]>) -> PyResult<Bound<'_, PyArray2<f32>>> {
    let rows = values.len();
    Array2::from_shape_vec((rows, 3), values.into_iter().flatten().collect())
        .map_err(value_error)
        .map(|value| value.into_pyarray(py))
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
