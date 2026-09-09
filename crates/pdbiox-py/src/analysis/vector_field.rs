//! Python vector-field values backed by the native integrator.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyclass(name = "StreamlineDirection", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyStreamlineDirection {
    Forward,
    Backward,
    Both,
}

impl From<PyStreamlineDirection> for pdbiox::analysis::StreamlineDirection {
    fn from(value: PyStreamlineDirection) -> Self {
        match value {
            PyStreamlineDirection::Forward => Self::Forward,
            PyStreamlineDirection::Backward => Self::Backward,
            PyStreamlineDirection::Both => Self::Both,
        }
    }
}

#[pyclass(name = "StreamlineOptions", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyStreamlineOptions {
    native: pdbiox::analysis::StreamlineOptions,
}

#[pymethods]
impl PyStreamlineOptions {
    #[new]
    #[pyo3(signature = (step_size=0.25, max_steps=4096, max_length=1000.0, min_speed=1.0e-6, direction=PyStreamlineDirection::Forward))]
    fn new(
        step_size: f32,
        max_steps: usize,
        max_length: f32,
        min_speed: f32,
        direction: PyStreamlineDirection,
    ) -> Self {
        Self {
            native: pdbiox::analysis::StreamlineOptions {
                step_size,
                max_steps,
                max_length,
                min_speed,
                direction: direction.into(),
            },
        }
    }

    #[getter]
    fn step_size(&self) -> f32 {
        self.native.step_size
    }
    #[getter]
    fn max_steps(&self) -> usize {
        self.native.max_steps
    }
    #[getter]
    fn max_length(&self) -> f32 {
        self.native.max_length
    }
    #[getter]
    fn min_speed(&self) -> f32 {
        self.native.min_speed
    }
}

#[pyclass(name = "VectorFieldGrid", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyVectorFieldGrid {
    native: pdbiox::analysis::VectorFieldGrid,
}

#[pymethods]
impl PyVectorFieldGrid {
    #[new]
    fn new(
        origin: [f32; 3],
        spacing: [f32; 3],
        shape: [usize; 3],
        vectors: Vec<[f32; 3]>,
    ) -> PyResult<Self> {
        let native = pdbiox::analysis::VectorFieldGrid::new(origin, spacing, shape, vectors)
            .map_err(value_error)?;
        Ok(Self { native })
    }

    #[getter]
    fn origin(&self) -> [f32; 3] {
        self.native.origin
    }
    #[getter]
    fn spacing(&self) -> [f32; 3] {
        self.native.spacing
    }
    #[getter]
    fn shape(&self) -> [usize; 3] {
        self.native.shape
    }
    fn sample(&self, point: [f32; 3]) -> Option<[f32; 3]> {
        self.native.sample(point)
    }
}

#[pyclass(name = "VectorFieldError", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyVectorFieldError {
    InvalidGrid,
    GridTooLarge,
    LengthMismatch,
    NonFiniteInput,
    InvalidOptions,
}

#[pyfunction]
fn integrate_streamlines(
    py: Python<'_>,
    field: &PyVectorFieldGrid,
    seeds: Vec<[f32; 3]>,
    options: &PyStreamlineOptions,
) -> PyResult<Vec<Vec<[f32; 3]>>> {
    py.detach(move || -> PyResult<Vec<Vec<[f32; 3]>>> {
        pdbiox::analysis::integrate_streamlines(&field.native, &seeds, options.native)
            .map_err(value_error)
    })
}

fn value_error(error: pdbiox::analysis::VectorFieldError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStreamlineDirection>()?;
    module.add_class::<PyStreamlineOptions>()?;
    module.add_class::<PyVectorFieldGrid>()?;
    module.add_class::<PyVectorFieldError>()?;
    module.add_function(wrap_pyfunction!(integrate_streamlines, module)?)
}
