//! Periodic, rotational and ordered-backbone geometry bindings.

use super::{PyEigenOptions, PyRotationMeanOptions};
use crate::intrinsic::PyRotation3;
use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "CircularSummary", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCircularSummary {
    #[pyo3(get)]
    mean: Option<f64>,
    #[pyo3(get)]
    resultant: f64,
    #[pyo3(get)]
    variance: f64,
}

#[pyclass(name = "BackboneFrame", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBackboneFrame {
    #[pyo3(get)]
    tangent: [f64; 3],
    #[pyo3(get)]
    normal: [f64; 3],
    #[pyo3(get)]
    binormal: [f64; 3],
    #[pyo3(get)]
    curvature: f64,
    #[pyo3(get)]
    torsion: Option<f64>,
}

#[pyclass(name = "HelixGeometry", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHelixGeometry {
    #[pyo3(get)]
    axis: [f64; 3],
    #[pyo3(get)]
    rise: f64,
    #[pyo3(get)]
    twist: f64,
}

#[pyclass(name = "BackboneTorsions", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBackboneTorsions {
    #[pyo3(get)]
    phi: Option<f64>,
    #[pyo3(get)]
    psi: Option<f64>,
    #[pyo3(get)]
    omega: Option<f64>,
}

#[pyclass(name = "BackboneCoordinates", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBackboneCoordinates {
    residues: Vec<molframe::BackboneResidue>,
}

#[pymethods]
impl PyBackboneCoordinates {
    #[new]
    fn new(
        py: Python<'_>,
        nitrogen: PyReadonlyArray2<'_, f32>,
        alpha_carbon: PyReadonlyArray2<'_, f32>,
        carbonyl_carbon: PyReadonlyArray2<'_, f32>,
        nitrogen_present: PyReadonlyArray1<'_, bool>,
        alpha_carbon_present: PyReadonlyArray1<'_, bool>,
        carbonyl_carbon_present: PyReadonlyArray1<'_, bool>,
        connected_to_next: PyReadonlyArray1<'_, bool>,
    ) -> PyResult<Self> {
        let (nitrogen, nitrogen_present) = optional_position_inputs(&nitrogen, &nitrogen_present)?;
        let (alpha_carbon, alpha_carbon_present) =
            optional_position_inputs(&alpha_carbon, &alpha_carbon_present)?;
        let (carbonyl_carbon, carbonyl_carbon_present) =
            optional_position_inputs(&carbonyl_carbon, &carbonyl_carbon_present)?;
        let connected = connected_to_next.as_slice().map_err(|_| {
            PyValueError::new_err(
                "connected_to_next must be C-contiguous; call numpy.ascontiguousarray explicitly",
            )
        })?;
        let count = nitrogen.len();
        if alpha_carbon.len() != count || carbonyl_carbon.len() != count || connected.len() != count
        {
            return Err(PyValueError::new_err("backbone arrays differ in length"));
        }
        let residues = py.detach(move || {
            (0..count)
                .map(|index| molframe::BackboneResidue {
                    nitrogen: nitrogen_present[index].then_some(nitrogen[index]),
                    alpha_carbon: alpha_carbon_present[index].then_some(alpha_carbon[index]),
                    carbonyl_carbon: carbonyl_carbon_present[index]
                        .then_some(carbonyl_carbon[index]),
                    connected_to_next: connected[index],
                })
                .collect()
        });
        Ok(Self { residues })
    }

    fn torsions(&self, py: Python<'_>) -> Vec<PyBackboneTorsions> {
        py.detach(|| molframe::backbone_torsions(&self.residues))
            .into_iter()
            .map(PyBackboneTorsions::from)
            .collect()
    }
}

#[pyfunction]
pub(crate) fn circular_summary(
    py: Python<'_>,
    radians: PyReadonlyArray1<'_, f64>,
    indeterminate_tolerance: f64,
) -> PyResult<PyCircularSummary> {
    let radians = radians.as_slice().map_err(|_| {
        PyValueError::new_err(
            "radians must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })?;
    py.detach(move || {
        let angles = periodic_angles(radians)?;
        molframe::circular_summary(&angles, indeterminate_tolerance)
    })
    .map(PyCircularSummary::from)
    .map_err(periodic_error)
}

#[pyfunction]
pub(crate) fn torus_summary(
    py: Python<'_>,
    radians: PyReadonlyArray2<'_, f64>,
    indeterminate_tolerance: f64,
) -> PyResult<Vec<PyCircularSummary>> {
    let shape = radians.shape();
    if shape.len() != 2 {
        return Err(PyValueError::new_err(
            "radians must have shape (n, dimensions)",
        ));
    }
    let dimensions = shape[1];
    if dimensions == 0 {
        return Err(PyValueError::new_err(
            "radians must have at least one angular dimension",
        ));
    }
    let radians = radians.as_slice().map_err(|_| {
        PyValueError::new_err(
            "radians must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })?;
    py.detach(move || {
        let points = radians
            .chunks_exact(dimensions)
            .map(|row| periodic_angles(row).map(Vec::into_boxed_slice))
            .collect::<Result<Vec<_>, _>>()?;
        molframe::torus_summary(&points, indeterminate_tolerance)
    })
    .map(|values| {
        values
            .iter()
            .copied()
            .map(PyCircularSummary::from)
            .collect()
    })
    .map_err(periodic_error)
}

#[pyfunction]
pub(crate) fn rotation_mean(
    py: Python<'_>,
    matrices: PyReadonlyArray3<'_, f64>,
    matrix_tolerance: f64,
    small_angle_tolerance: f64,
    half_turn_tolerance: f64,
    mean_tolerance: f64,
    maximum_iterations: usize,
) -> PyResult<PyRotation3> {
    let shape = matrices.shape();
    if shape.len() != 3 || shape[1] != 3 || shape[2] != 3 {
        return Err(PyValueError::new_err("matrices must have shape (n, 3, 3)"));
    }
    let matrices = matrix_values(&matrices)?;
    let options = molframe::RotationMeanOptions {
        convergence_tolerance: mean_tolerance,
        maximum_iterations,
        rotation: molframe::RotationOptions {
            matrix_tolerance,
            small_angle_tolerance,
            half_turn_tolerance,
        },
    };
    py.detach(move || {
        let rotations = rotations_from_matrices(matrices, options.rotation.matrix_tolerance)?;
        molframe::rotation_mean_with_options(&rotations, options)
    })
    .map(PyRotation3::from_rust)
    .map_err(rotation_error)
}

#[pyfunction]
pub(crate) fn rotation_mean_with_options(
    py: Python<'_>,
    matrices: PyReadonlyArray3<'_, f64>,
    options: &PyRotationMeanOptions,
) -> PyResult<PyRotation3> {
    let shape = matrices.shape();
    if shape.len() != 3 || shape[1] != 3 || shape[2] != 3 {
        return Err(PyValueError::new_err("matrices must have shape (n, 3, 3)"));
    }
    let matrices = matrix_values(&matrices)?;
    let options = options.0;
    py.detach(move || {
        let rotations = rotations_from_matrices(matrices, options.rotation.matrix_tolerance)?;
        molframe::rotation_mean_with_options(&rotations, options)
    })
    .map(PyRotation3::from_rust)
    .map_err(rotation_error)
}

#[pyfunction]
pub(crate) fn path_torsions(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    present: PyReadonlyArray1<'_, bool>,
    maximum: usize,
) -> PyResult<Vec<Option<f64>>> {
    let (positions, present) = optional_position_inputs(&positions, &present)?;
    Ok(py.detach(move || {
        let positions = materialize_optional_positions(positions, present);
        molframe::path_torsions(&positions, maximum)
    }))
}

#[pyfunction]
pub(crate) fn backbone_frames(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    present: PyReadonlyArray1<'_, bool>,
) -> PyResult<Vec<Option<PyBackboneFrame>>> {
    let (positions, present) = optional_position_inputs(&positions, &present)?;
    Ok(py
        .detach(move || {
            let positions = materialize_optional_positions(positions, present);
            molframe::backbone_frames(&positions)
        })
        .into_iter()
        .map(|value| value.map(Into::into))
        .collect())
}

#[pyfunction]
pub(crate) fn helix_geometry(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    present: PyReadonlyArray1<'_, bool>,
    options: &PyEigenOptions,
) -> PyResult<Option<PyHelixGeometry>> {
    let (positions, present) = optional_position_inputs(&positions, &present)?;
    let options = options.inner;
    py.detach(move || {
        let positions = materialize_optional_positions(positions, present);
        molframe::helix_geometry_with_options(&positions, options)
    })
    .map(|value| value.map(PyHelixGeometry::from))
    .map_err(eigen_error)
}

#[pyfunction]
pub(crate) fn helix_geometry_with_options(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    present: PyReadonlyArray1<'_, bool>,
    options: &PyEigenOptions,
) -> PyResult<Option<PyHelixGeometry>> {
    let (positions, present) = optional_position_inputs(&positions, &present)?;
    let options = options.inner;
    py.detach(move || {
        let positions = materialize_optional_positions(positions, present);
        molframe::helix_geometry_with_options(&positions, options)
    })
    .map(|value| value.map(PyHelixGeometry::from))
    .map_err(eigen_error)
}

fn periodic_angles(
    values: &[f64],
) -> Result<Vec<molframe::PeriodicAngle>, molframe::PeriodicError> {
    values
        .iter()
        .copied()
        .map(molframe::PeriodicAngle::from_radians)
        .collect()
}

fn optional_position_inputs<'a>(
    positions: &'a PyReadonlyArray2<'_, f32>,
    present_array: &'a PyReadonlyArray1<'_, bool>,
) -> PyResult<(&'a [[f32; 3]], &'a [bool])> {
    if positions.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("positions must have shape (n, 3)"));
    }
    let values = positions.as_slice().map_err(|_| {
        PyValueError::new_err(
            "positions must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })?;
    let (positions, remainder) = values.as_chunks::<3>();
    if !remainder.is_empty() {
        return Err(PyValueError::new_err("positions must have shape (n, 3)"));
    }
    let present = present_array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "present must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })?;
    if positions.len() != present.len() {
        return Err(PyValueError::new_err(
            "positions and present differ in length",
        ));
    }
    Ok((positions, present))
}

fn materialize_optional_positions(
    positions: &[[f32; 3]],
    present: &[bool],
) -> Vec<Option<[f32; 3]>> {
    positions
        .iter()
        .zip(present)
        .map(|(&position, &present)| present.then_some(position))
        .collect()
}

fn matrix_values<'a>(matrices: &'a PyReadonlyArray3<'_, f64>) -> PyResult<&'a [f64]> {
    matrices.as_slice().map_err(|_| {
        PyValueError::new_err(
            "matrices must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })
}

fn rotations_from_matrices(
    matrices: &[f64],
    matrix_tolerance: f64,
) -> Result<Vec<molframe::Rotation3>, molframe::RotationError> {
    matrices
        .as_chunks::<9>()
        .0
        .iter()
        .map(|matrix| {
            molframe::Rotation3::from_matrix(
                [
                    [matrix[0], matrix[1], matrix[2]],
                    [matrix[3], matrix[4], matrix[5]],
                    [matrix[6], matrix[7], matrix[8]],
                ],
                matrix_tolerance,
            )
        })
        .collect()
}

fn periodic_error(error: molframe::PeriodicError) -> PyErr {
    PyValueError::new_err(format!("invalid periodic data: {error:?}"))
}

fn rotation_error(error: molframe::RotationError) -> PyErr {
    PyValueError::new_err(format!("invalid rotation data: {error:?}"))
}

fn eigen_error(error: molframe::EigenError) -> PyErr {
    PyValueError::new_err(format!("eigendecomposition failed: {error:?}"))
}

impl From<molframe::CircularSummary> for PyCircularSummary {
    fn from(value: molframe::CircularSummary) -> Self {
        Self {
            mean: value.mean.map(molframe::PeriodicAngle::radians),
            resultant: value.resultant,
            variance: value.variance,
        }
    }
}

impl From<molframe::BackboneFrame> for PyBackboneFrame {
    fn from(value: molframe::BackboneFrame) -> Self {
        Self {
            tangent: value.tangent,
            normal: value.normal,
            binormal: value.binormal,
            curvature: value.curvature,
            torsion: value.torsion,
        }
    }
}

impl From<molframe::HelixGeometry> for PyHelixGeometry {
    fn from(value: molframe::HelixGeometry) -> Self {
        Self {
            axis: value.axis,
            rise: value.rise,
            twist: value.twist,
        }
    }
}

impl From<molframe::BackboneTorsions> for PyBackboneTorsions {
    fn from(value: molframe::BackboneTorsions) -> Self {
        Self {
            phi: value.phi,
            psi: value.psi,
            omega: value.omega,
        }
    }
}
