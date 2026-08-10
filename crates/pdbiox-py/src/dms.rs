//! Typed Python view of the native DMS topology-and-frame model.

use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pdbiox::traj::{DmsParticle, DmsSystem, read_dms, write_dms};
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;

const VECTOR_WIDTH: usize = 3;
const BOND_WIDTH: usize = 2;
const MASK_ARGUMENT: &str = "mask";
const NUMPY_MASKED_ARRAY_MODULE: &str = "numpy.ma";
const ARRAY_FUNCTION: &str = "array";
const SET_FLAGS_METHOD: &str = "setflags";
const WRITE_ARGUMENT: &str = "write";

#[pyclass(name = "DmsParticle", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyDmsParticle {
    value: DmsParticle,
}

#[pymethods]
impl PyDmsParticle {
    #[getter]
    const fn atomic_number(&self) -> Option<u16> {
        self.value.atomic_number
    }
    #[getter]
    const fn component(&self) -> Option<i64> {
        self.value.component
    }
    #[getter]
    const fn nonbonded_type(&self) -> Option<i64> {
        self.value.nonbonded_type
    }
    #[getter]
    const fn mass(&self) -> Option<f64> {
        self.value.mass
    }
    #[getter]
    const fn charge(&self) -> Option<f64> {
        self.value.charge
    }
    #[getter]
    const fn residue_id(&self) -> Option<i64> {
        self.value.residue_id
    }
    #[getter]
    fn residue_name(&self) -> Option<&str> {
        self.value.residue_name.as_deref()
    }
    #[getter]
    fn chain(&self) -> Option<&str> {
        self.value.chain.as_deref()
    }
    #[getter]
    fn segment(&self) -> Option<&str> {
        self.value.segment.as_deref()
    }
    #[getter]
    fn name(&self) -> Option<&str> {
        self.value.name.as_deref()
    }
    #[getter]
    fn insertion(&self) -> Option<&str> {
        self.value.insertion.as_deref()
    }
    #[getter]
    const fn formal_charge(&self) -> Option<f64> {
        self.value.formal_charge
    }
    #[getter]
    const fn occupancy(&self) -> Option<f64> {
        self.value.occupancy
    }
    #[getter]
    const fn b_factor(&self) -> Option<f64> {
        self.value.b_factor
    }
    #[getter]
    const fn temperature_group(&self) -> Option<i64> {
        self.value.temperature_group
    }
    #[getter]
    const fn energy_group(&self) -> Option<i64> {
        self.value.energy_group
    }
    #[getter]
    const fn ligand_group(&self) -> Option<i64> {
        self.value.ligand_group
    }
    #[getter]
    const fn bias_group(&self) -> Option<i64> {
        self.value.bias_group
    }
}

#[pyclass(name = "DmsSystem", frozen, skip_from_py_object)]
pub(crate) struct PyDmsSystem {
    inner: DmsSystem,
    positions: Py<PyArray2<f64>>,
    velocities: Py<PyAny>,
    bond_indices: Py<PyArray2<u32>>,
    bond_orders: Py<PyArray1<f64>>,
    cell: Option<Py<PyArray2<f64>>>,
}

#[pymethods]
impl PyDmsSystem {
    #[staticmethod]
    fn read(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        read_dms(path)
            .map_err(dms_error)
            .and_then(|system| Self::from_system(py, system))
    }

    fn write(&self, path: PathBuf) -> PyResult<()> {
        write_dms(path, &self.inner).map_err(dms_error)
    }

    fn __len__(&self) -> usize {
        self.inner.topology.particles.len()
    }

    #[getter]
    fn version(&self) -> Option<(u32, u32)> {
        self.inner
            .version
            .map(|version| (version.major, version.minor))
    }

    #[getter]
    fn particles(&self) -> Vec<PyDmsParticle> {
        self.inner
            .topology
            .particles
            .iter()
            .cloned()
            .map(|value| PyDmsParticle { value })
            .collect()
    }

    #[getter]
    fn positions(&self, py: Python<'_>) -> Py<PyArray2<f64>> {
        self.positions.clone_ref(py)
    }

    #[getter]
    fn velocities(&self, py: Python<'_>) -> Py<PyAny> {
        self.velocities.clone_ref(py)
    }

    #[getter]
    fn bond_indices(&self, py: Python<'_>) -> Py<PyArray2<u32>> {
        self.bond_indices.clone_ref(py)
    }

    #[getter]
    fn bond_orders(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.bond_orders.clone_ref(py)
    }

    #[getter]
    fn cell(&self, py: Python<'_>) -> Option<Py<PyArray2<f64>>> {
        self.cell.as_ref().map(|cell| cell.clone_ref(py))
    }
}

impl PyDmsSystem {
    fn from_system(py: Python<'_>, inner: DmsSystem) -> PyResult<Self> {
        let positions = vector_array(py, &inner.frame.positions)?;
        let velocities = masked_velocities(py, &inner.frame.velocities)?;
        let bond_indices = Array2::from_shape_vec(
            (inner.topology.bonds.len(), BOND_WIDTH),
            inner
                .topology
                .bonds
                .iter()
                .flat_map(|bond| [bond.p0, bond.p1])
                .collect(),
        )
        .map_err(|error| shape_error(&error))?
        .into_pyarray(py);
        let _readonly = bond_indices.readwrite().make_nonwriteable();
        let bond_indices = bond_indices.unbind();
        let bond_orders =
            Array1::from_vec(inner.topology.bonds.iter().map(|bond| bond.order).collect())
                .into_pyarray(py);
        let _readonly = bond_orders.readwrite().make_nonwriteable();
        let bond_orders = bond_orders.unbind();
        let cell = inner
            .frame
            .cell
            .map(|cell| vector_array(py, &cell.vectors))
            .transpose()?;
        Ok(Self {
            inner,
            positions,
            velocities,
            bond_indices,
            bond_orders,
            cell,
        })
    }
}

fn vector_array<const N: usize>(
    py: Python<'_>,
    values: &[[f64; N]],
) -> PyResult<Py<PyArray2<f64>>> {
    Array2::from_shape_vec(
        (values.len(), N),
        values
            .iter()
            .flat_map(|value| value.iter().copied())
            .collect(),
    )
    .map_err(|error| shape_error(&error))
    .map(|array| {
        let array = array.into_pyarray(py);
        let _readonly = array.readwrite().make_nonwriteable();
        array.unbind()
    })
}

fn masked_velocities(
    py: Python<'_>,
    values: &[Option<[f64; VECTOR_WIDTH]>],
) -> PyResult<Py<PyAny>> {
    let mut vectors = Vec::with_capacity(values.len() * VECTOR_WIDTH);
    let mut mask = Vec::with_capacity(values.len() * VECTOR_WIDTH);
    for value in values {
        if let Some(value) = value {
            vectors.extend(value);
            mask.extend([false; VECTOR_WIDTH]);
        } else {
            vectors.extend([0.0; VECTOR_WIDTH]);
            mask.extend([true; VECTOR_WIDTH]);
        }
    }
    let vectors = Array2::from_shape_vec((values.len(), VECTOR_WIDTH), vectors)
        .map_err(|error| shape_error(&error))?
        .into_pyarray(py);
    let mask = Array2::from_shape_vec((values.len(), VECTOR_WIDTH), mask)
        .map_err(|error| shape_error(&error))?
        .into_pyarray(py);
    let kwargs = PyDict::new(py);
    kwargs.set_item(MASK_ARGUMENT, mask)?;
    let result = py
        .import(NUMPY_MASKED_ARRAY_MODULE)?
        .getattr(ARRAY_FUNCTION)?
        .call((vectors,), Some(&kwargs))?;
    let flags = PyDict::new(py);
    flags.set_item(WRITE_ARGUMENT, false)?;
    result.call_method(SET_FLAGS_METHOD, (), Some(&flags))?;
    Ok(result.unbind())
}

fn dms_error(error: pdbiox::traj::DmsError) -> PyErr {
    match error {
        pdbiox::traj::DmsError::Io(error) => PyOSError::new_err(error.to_string()),
        error => PyValueError::new_err(error.to_string()),
    }
}

fn shape_error(error: &numpy::ndarray::ShapeError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[cfg(test)]
#[path = "dms_tests.rs"]
mod tests;
