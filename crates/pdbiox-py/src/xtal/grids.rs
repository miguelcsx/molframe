//! Native text-grid and bounded MRC block-reader projections.

use std::path::PathBuf;

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::crystallography::PyUnitCell;
use super::maps::PyDensityMap;

#[pyclass(name = "CubeAtom", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCubeAtom {
    #[pyo3(get)]
    number: i32,
    #[pyo3(get)]
    charge: f64,
    #[pyo3(get)]
    position: [f64; 3],
}

impl From<pdbiox::xtal::CubeAtom> for PyCubeAtom {
    fn from(value: pdbiox::xtal::CubeAtom) -> Self {
        Self {
            number: value.number,
            charge: value.charge,
            position: value.position,
        }
    }
}

#[pyclass(name = "CubeGrid", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCubeGrid {
    map: PyDensityMap,
    atoms: Vec<PyCubeAtom>,
    orbitals: Vec<i32>,
    fields: usize,
}

#[pymethods]
impl PyCubeGrid {
    #[getter]
    fn map(&self) -> PyDensityMap {
        self.map.clone()
    }
    #[getter]
    fn atoms(&self) -> Vec<PyCubeAtom> {
        self.atoms.clone()
    }
    #[getter]
    fn orbitals(&self) -> &[i32] {
        &self.orbitals
    }
    #[getter]
    fn fields(&self) -> usize {
        self.fields
    }
}

impl From<pdbiox::xtal::CubeGrid> for PyCubeGrid {
    fn from(value: pdbiox::xtal::CubeGrid) -> Self {
        Self {
            map: PyDensityMap(value.map),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            orbitals: value.orbitals,
            fields: value.fields,
        }
    }
}

#[pyclass(name = "MrcBlockOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMrcBlockOptions(pdbiox::xtal::MrcBlockOptions);

#[pymethods]
impl PyMrcBlockOptions {
    #[new]
    #[pyo3(signature = (*, memory_limit_bytes=pdbiox::xtal::DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES))]
    fn new(memory_limit_bytes: usize) -> Self {
        Self(pdbiox::xtal::MrcBlockOptions::default().with_memory_limit(memory_limit_bytes))
    }
    #[getter]
    fn memory_limit_bytes(&self) -> usize {
        self.0.memory_limit_bytes
    }
}

#[pyclass(name = "MrcMapDescriptor", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMrcMapDescriptor(pdbiox::xtal::MrcMapDescriptor);

#[pymethods]
impl PyMrcMapDescriptor {
    #[getter]
    fn dimensions(&self) -> [usize; 3] {
        self.0.dimensions
    }
    #[getter]
    fn starts(&self) -> [i32; 3] {
        self.0.starts
    }
    #[getter]
    fn sampling(&self) -> [usize; 3] {
        self.0.sampling
    }
    #[getter]
    fn cell(&self) -> PyResult<PyUnitCell> {
        PyUnitCell::from_native(self.0.cell)
    }
    #[getter]
    fn origin(&self) -> [f64; 3] {
        self.0.origin
    }
    #[getter]
    fn space_group(&self) -> i32 {
        self.0.space_group
    }
    #[getter]
    fn labels(&self) -> Vec<&str> {
        self.0.labels.iter().map(AsRef::as_ref).collect()
    }
    #[getter]
    fn extended_header_bytes(&self) -> usize {
        self.0.extended_header_bytes
    }
    #[getter]
    fn mode(&self) -> i32 {
        self.0.mode
    }
    #[getter]
    fn value_range(&self) -> [f32; 2] {
        self.0.value_range
    }
}

#[pyclass(name = "MrcBlockReader", skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyMrcBlockReader {
    inner: pdbiox::xtal::MrcBlockReader<std::fs::File>,
    values: Vec<f32>,
}

#[pymethods]
impl PyMrcBlockReader {
    #[new]
    #[pyo3(signature = (path, *, options=None))]
    fn new(path: PathBuf, options: Option<PyMrcBlockOptions>) -> PyResult<Self> {
        let options = match options {
            Some(value) => value.0,
            None => pdbiox::xtal::MrcBlockOptions::default(),
        };
        pdbiox::xtal::MrcBlockReader::open(path, options)
            .map(|inner| Self {
                inner,
                values: Vec::new(),
            })
            .map_err(value_error)
    }

    #[getter]
    fn descriptor(&self) -> PyMrcMapDescriptor {
        PyMrcMapDescriptor(self.inner.descriptor().clone())
    }

    fn read_block<'py>(
        &mut self,
        py: Python<'py>,
        origin: [usize; 3],
        dimensions: [usize; 3],
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        self.inner
            .read_block_into(origin, dimensions, &mut self.values)
            .map_err(value_error)?;
        Ok(Array1::from_vec(self.values.clone()).into_pyarray(py))
    }
}

#[pyfunction]
pub(crate) fn read_cube(py: Python<'_>, data: &[u8]) -> PyResult<PyCubeGrid> {
    py.detach(move || -> PyResult<PyCubeGrid> {
        pdbiox::xtal::read_cube(data)
            .map(Into::into)
            .map_err(value_error)
    })
}

#[pyfunction]
pub(crate) fn read_dx(py: Python<'_>, data: &[u8]) -> PyResult<PyDensityMap> {
    py.detach(move || -> PyResult<PyDensityMap> {
        pdbiox::xtal::read_dx(data)
            .map(PyDensityMap)
            .map_err(value_error)
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES",
        pdbiox::xtal::DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES,
    )?;
    module.add_class::<PyCubeAtom>()?;
    module.add_class::<PyCubeGrid>()?;
    module.add_class::<PyMrcBlockOptions>()?;
    module.add_class::<PyMrcMapDescriptor>()?;
    module.add_class::<PyMrcBlockReader>()?;
    module.add_function(wrap_pyfunction!(read_cube, module)?)?;
    module.add_function(wrap_pyfunction!(read_dx, module)?)?;
    Ok(())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
