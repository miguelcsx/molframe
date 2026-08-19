//! Typed adapters for public geometry controls and small native results.

use super::moments::PyEigenOptions;
use ::pdbiox;
use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, EigenError, PyValueError);
create_exception!(_native, FluctuationError, PyValueError);
create_exception!(_native, MatrixError, PyValueError);
create_exception!(_native, PeriodicError, PyValueError);
create_exception!(_native, RotationError, PyValueError);
create_exception!(_native, SuperposeError, PyValueError);

#[derive(Clone, Copy, Debug)]
enum DecompositionValue {
    Three(pdbiox::Decomposition<3>),
    Four(pdbiox::Decomposition<4>),
}

#[pyclass(name = "Decomposition", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDecomposition {
    value: DecompositionValue,
}

#[pymethods]
impl PyDecomposition {
    #[getter]
    fn values(&self) -> Vec<f64> {
        match self.value {
            DecompositionValue::Three(value) => value.values.to_vec(),
            DecompositionValue::Four(value) => value.values.to_vec(),
        }
    }

    #[getter]
    fn vectors(&self) -> Vec<Vec<f64>> {
        match self.value {
            DecompositionValue::Three(value) => {
                value.vectors.iter().map(|row| row.to_vec()).collect()
            }
            DecompositionValue::Four(value) => {
                value.vectors.iter().map(|row| row.to_vec()).collect()
            }
        }
    }

    fn dominant(&self) -> Vec<f64> {
        match self.value {
            DecompositionValue::Three(value) => value.dominant().to_vec(),
            DecompositionValue::Four(value) => value.dominant().to_vec(),
        }
    }

    fn vector(&self, position: usize) -> Option<Vec<f64>> {
        match self.value {
            DecompositionValue::Three(value) => {
                value.vector(position).map(|vector| vector.to_vec())
            }
            DecompositionValue::Four(value) => value.vector(position).map(|vector| vector.to_vec()),
        }
    }
}

#[pyclass(name = "TorusMetric", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTorusMetric(pub(crate) pdbiox::TorusMetric);

#[pymethods]
impl PyTorusMetric {
    #[new]
    fn new(weights: Vec<f64>) -> PyResult<Self> {
        pdbiox::TorusMetric::new(weights.into_boxed_slice())
            .map(Self)
            .map_err(periodic_error)
    }

    #[staticmethod]
    fn uniform(dimension: usize) -> PyResult<Self> {
        pdbiox::TorusMetric::uniform(dimension)
            .map(Self)
            .map_err(periodic_error)
    }

    #[getter]
    fn dimension(&self) -> usize {
        self.0.dimension()
    }

    fn distance(
        &self,
        left: PyReadonlyArray1<'_, f64>,
        right: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<f64> {
        let left = periodic_values(left.as_slice()?)?;
        let right = periodic_values(right.as_slice()?)?;
        self.0.distance(&left, &right).map_err(periodic_error)
    }
}

#[pyclass(name = "BackboneResidue", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBackboneResidue(pub(crate) pdbiox::BackboneResidue);

#[pymethods]
impl PyBackboneResidue {
    #[new]
    #[pyo3(signature = (nitrogen=None, alpha_carbon=None, carbonyl_carbon=None, connected_to_next=false))]
    fn new(
        nitrogen: Option<[f32; 3]>,
        alpha_carbon: Option<[f32; 3]>,
        carbonyl_carbon: Option<[f32; 3]>,
        connected_to_next: bool,
    ) -> PyResult<Self> {
        if [nitrogen, alpha_carbon, carbonyl_carbon]
            .into_iter()
            .flatten()
            .flatten()
            .any(|value| !value.is_finite())
        {
            return Err(PyValueError::new_err("backbone coordinates must be finite"));
        }
        Ok(Self(pdbiox::BackboneResidue {
            nitrogen,
            alpha_carbon,
            carbonyl_carbon,
            connected_to_next,
        }))
    }

    #[getter]
    const fn nitrogen(&self) -> Option<[f32; 3]> {
        self.0.nitrogen
    }

    #[getter]
    const fn alpha_carbon(&self) -> Option<[f32; 3]> {
        self.0.alpha_carbon
    }

    #[getter]
    const fn carbonyl_carbon(&self) -> Option<[f32; 3]> {
        self.0.carbonyl_carbon
    }

    #[getter]
    const fn connected_to_next(&self) -> bool {
        self.0.connected_to_next
    }
}

#[pyclass(name = "RotationOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyRotationOptions(pub(crate) pdbiox::RotationOptions);

#[pymethods]
impl PyRotationOptions {
    #[classattr]
    const STANDARD_MATRIX_TOLERANCE: f64 = pdbiox::RotationOptions::STANDARD_MATRIX_TOLERANCE;

    #[classattr]
    const STANDARD_SMALL_ANGLE_TOLERANCE: f64 =
        pdbiox::RotationOptions::STANDARD_SMALL_ANGLE_TOLERANCE;

    #[classattr]
    const STANDARD_HALF_TURN_TOLERANCE: f64 = pdbiox::RotationOptions::STANDARD_HALF_TURN_TOLERANCE;

    #[new]
    fn new(matrix_tolerance: f64, small_angle_tolerance: f64, half_turn_tolerance: f64) -> Self {
        Self(pdbiox::RotationOptions {
            matrix_tolerance,
            small_angle_tolerance,
            half_turn_tolerance,
        })
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(pdbiox::RotationOptions::standard())
    }

    #[getter]
    const fn matrix_tolerance(&self) -> f64 {
        self.0.matrix_tolerance
    }

    #[getter]
    const fn small_angle_tolerance(&self) -> f64 {
        self.0.small_angle_tolerance
    }

    #[getter]
    const fn half_turn_tolerance(&self) -> f64 {
        self.0.half_turn_tolerance
    }
}

#[pyclass(name = "RotationMeanOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyRotationMeanOptions(pub(crate) pdbiox::RotationMeanOptions);

#[pymethods]
impl PyRotationMeanOptions {
    #[new]
    fn new(
        convergence_tolerance: f64,
        maximum_iterations: usize,
        rotation: &PyRotationOptions,
    ) -> Self {
        Self(pdbiox::RotationMeanOptions {
            convergence_tolerance,
            maximum_iterations,
            rotation: rotation.0,
        })
    }

    #[getter]
    const fn convergence_tolerance(&self) -> f64 {
        self.0.convergence_tolerance
    }

    #[getter]
    const fn maximum_iterations(&self) -> usize {
        self.0.maximum_iterations
    }

    #[getter]
    const fn rotation(&self) -> PyRotationOptions {
        PyRotationOptions(self.0.rotation)
    }
}

#[pyclass(name = "SuperposeOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySuperposeOptions(pub(crate) pdbiox::SuperposeOptions);

#[pymethods]
impl PySuperposeOptions {
    #[classattr]
    const STANDARD_COLLINEAR_RELATIVE_TOLERANCE: f64 =
        pdbiox::SuperposeOptions::STANDARD_COLLINEAR_RELATIVE_TOLERANCE;

    #[new]
    fn new(collinear_relative_tolerance: f64, eigen: &PyEigenOptions) -> Self {
        Self(pdbiox::SuperposeOptions {
            collinear_relative_tolerance,
            eigen: eigen.inner,
        })
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(pdbiox::SuperposeOptions::standard())
    }

    #[getter]
    const fn collinear_relative_tolerance(&self) -> f64 {
        self.0.collinear_relative_tolerance
    }

    #[getter]
    const fn eigen(&self) -> PyEigenOptions {
        PyEigenOptions {
            inner: self.0.eigen,
        }
    }
}

#[pyfunction]
pub(crate) fn symmetric(
    py: Python<'_>,
    matrix: PyReadonlyArray2<'_, f64>,
) -> PyResult<PyDecomposition> {
    symmetric_impl(py, matrix, None)
}

#[pyfunction]
pub(crate) fn symmetric_with_options(
    py: Python<'_>,
    matrix: PyReadonlyArray2<'_, f64>,
    options: &PyEigenOptions,
) -> PyResult<PyDecomposition> {
    symmetric_impl(py, matrix, Some(options.inner))
}

fn symmetric_impl(
    py: Python<'_>,
    matrix: PyReadonlyArray2<'_, f64>,
    options: Option<pdbiox::EigenOptions>,
) -> PyResult<PyDecomposition> {
    let values = square_values(matrix)?;
    match values.len() {
        9 => {
            let matrix = array3(&values);
            py.detach(move || match options {
                Some(options) => pdbiox::geom::symmetric_with_options(matrix, options),
                None => pdbiox::geom::symmetric(matrix),
            })
            .map(|value| PyDecomposition {
                value: DecompositionValue::Three(value),
            })
            .map_err(eigen_error)
        }
        16 => {
            let matrix = array4(&values);
            py.detach(move || match options {
                Some(options) => pdbiox::geom::symmetric_with_options(matrix, options),
                None => pdbiox::geom::symmetric(matrix),
            })
            .map(|value| PyDecomposition {
                value: DecompositionValue::Four(value),
            })
            .map_err(eigen_error)
        }
        _ => Err(PyValueError::new_err(
            "symmetric decomposition supports only 3x3 or 4x4 matrices",
        )),
    }
}

#[pyfunction]
pub(crate) fn backbone_torsions(
    py: Python<'_>,
    residues: Vec<PyBackboneResidue>,
) -> Vec<super::intrinsic::PyBackboneTorsions> {
    py.detach(move || {
        pdbiox::backbone_torsions(&residues.iter().map(|value| value.0).collect::<Vec<_>>())
    })
    .into_iter()
    .map(super::intrinsic::PyBackboneTorsions::from)
    .collect()
}

fn square_values(matrix: PyReadonlyArray2<'_, f64>) -> PyResult<Vec<f64>> {
    let shape = matrix.shape();
    if shape.len() != 2 || shape[0] != shape[1] || !matches!(shape[0], 3 | 4) {
        return Err(PyValueError::new_err(
            "matrix must have shape (3, 3) or (4, 4)",
        ));
    }
    let values = matrix.as_slice()?.to_vec();
    drop(matrix);
    Ok(values)
}

fn array3(values: &[f64]) -> [[f64; 3]; 3] {
    [
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]],
    ]
}

fn array4(values: &[f64]) -> [[f64; 4]; 4] {
    [
        [values[0], values[1], values[2], values[3]],
        [values[4], values[5], values[6], values[7]],
        [values[8], values[9], values[10], values[11]],
        [values[12], values[13], values[14], values[15]],
    ]
}

fn periodic_values(values: &[f64]) -> PyResult<Vec<pdbiox::PeriodicAngle>> {
    values
        .iter()
        .copied()
        .map(|value| pdbiox::PeriodicAngle::from_radians(value).map_err(periodic_error))
        .collect()
}

fn eigen_error(error: pdbiox::EigenError) -> PyErr {
    EigenError::new_err(format!("{error:?}"))
}

fn periodic_error(error: pdbiox::PeriodicError) -> PyErr {
    PeriodicError::new_err(format!("{error:?}"))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    for (name, value) in [
        ("EigenError", module.py().get_type::<EigenError>()),
        (
            "FluctuationError",
            module.py().get_type::<FluctuationError>(),
        ),
        ("MatrixError", module.py().get_type::<MatrixError>()),
        ("PeriodicError", module.py().get_type::<PeriodicError>()),
        ("RotationError", module.py().get_type::<RotationError>()),
        ("SuperposeError", module.py().get_type::<SuperposeError>()),
    ] {
        module.add(name, value)?;
    }
    module.add_class::<PyDecomposition>()?;
    super::distance_types::register(module)?;
    module.add_class::<PyTorusMetric>()?;
    module.add_class::<PyBackboneResidue>()?;
    module.add_class::<PyRotationOptions>()?;
    module.add_class::<PyRotationMeanOptions>()?;
    module.add_class::<PySuperposeOptions>()?;
    Ok(())
}
