//! Periodic, rotational and ordered-backbone geometry bindings.

use super::PyEigenOptions;
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
    residues: Vec<pdbiox::BackboneResidue>,
}

#[pymethods]
impl PyBackboneCoordinates {
    #[new]
    fn new(
        nitrogen: PyReadonlyArray2<'_, f32>,
        alpha_carbon: PyReadonlyArray2<'_, f32>,
        carbonyl_carbon: PyReadonlyArray2<'_, f32>,
        nitrogen_present: PyReadonlyArray1<'_, bool>,
        alpha_carbon_present: PyReadonlyArray1<'_, bool>,
        carbonyl_carbon_present: PyReadonlyArray1<'_, bool>,
        connected_to_next: PyReadonlyArray1<'_, bool>,
    ) -> PyResult<Self> {
        let nitrogen = optional_positions(nitrogen, nitrogen_present)?;
        let alpha_carbon = optional_positions(alpha_carbon, alpha_carbon_present)?;
        let carbonyl_carbon = optional_positions(carbonyl_carbon, carbonyl_carbon_present)?;
        let connected = connected_to_next.as_slice()?.to_vec();
        drop(connected_to_next);
        let count = nitrogen.len();
        if alpha_carbon.len() != count || carbonyl_carbon.len() != count || connected.len() != count
        {
            return Err(PyValueError::new_err("backbone arrays differ in length"));
        }
        let residues = (0..count)
            .map(|index| pdbiox::BackboneResidue {
                nitrogen: nitrogen[index],
                alpha_carbon: alpha_carbon[index],
                carbonyl_carbon: carbonyl_carbon[index],
                connected_to_next: connected[index],
            })
            .collect();
        Ok(Self { residues })
    }

    fn torsions(&self, py: Python<'_>) -> Vec<PyBackboneTorsions> {
        py.detach(|| pdbiox::backbone_torsions(&self.residues))
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
    let angles = periodic_angles(radians.as_slice()?)?;
    drop(radians);
    py.detach(|| pdbiox::circular_summary(&angles, indeterminate_tolerance))
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
    let points = radians
        .as_array()
        .outer_iter()
        .map(|row| {
            let values = row.iter().copied().collect::<Vec<_>>();
            periodic_angles(&values).map(Vec::into_boxed_slice)
        })
        .collect::<PyResult<Vec<_>>>()?;
    drop(radians);
    py.detach(|| pdbiox::torus_summary(&points, indeterminate_tolerance))
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
    let rotations = matrices
        .as_array()
        .outer_iter()
        .map(|matrix| {
            pdbiox::Rotation3::from_matrix(
                [
                    [matrix[[0, 0]], matrix[[0, 1]], matrix[[0, 2]]],
                    [matrix[[1, 0]], matrix[[1, 1]], matrix[[1, 2]]],
                    [matrix[[2, 0]], matrix[[2, 1]], matrix[[2, 2]]],
                ],
                matrix_tolerance,
            )
            .map_err(rotation_error)
        })
        .collect::<PyResult<Vec<_>>>()?;
    drop(matrices);
    let options = pdbiox::RotationMeanOptions {
        convergence_tolerance: mean_tolerance,
        maximum_iterations,
        rotation: pdbiox::RotationOptions {
            matrix_tolerance,
            small_angle_tolerance,
            half_turn_tolerance,
        },
    };
    py.detach(|| pdbiox::rotation_mean_with_options(&rotations, options))
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
    let positions = optional_positions(positions, present)?;
    Ok(py.detach(|| pdbiox::path_torsions(&positions, maximum)))
}

#[pyfunction]
pub(crate) fn backbone_frames(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    present: PyReadonlyArray1<'_, bool>,
) -> PyResult<Vec<Option<PyBackboneFrame>>> {
    let positions = optional_positions(positions, present)?;
    Ok(py
        .detach(|| pdbiox::backbone_frames(&positions))
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
    let positions = optional_positions(positions, present)?;
    py.detach(|| pdbiox::helix_geometry_with_options(&positions, options.inner))
        .map(|value| value.map(PyHelixGeometry::from))
        .map_err(eigen_error)
}

fn periodic_angles(values: &[f64]) -> PyResult<Vec<pdbiox::PeriodicAngle>> {
    values
        .iter()
        .copied()
        .map(|value| pdbiox::PeriodicAngle::from_radians(value).map_err(periodic_error))
        .collect()
}

fn optional_positions(
    positions: PyReadonlyArray2<'_, f32>,
    present_array: PyReadonlyArray1<'_, bool>,
) -> PyResult<Vec<Option<[f32; 3]>>> {
    if positions.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("positions must have shape (n, 3)"));
    }
    let present = present_array.as_slice()?.to_vec();
    drop(present_array);
    let rows = positions.as_array();
    if rows.nrows() != present.len() {
        return Err(PyValueError::new_err(
            "positions and present differ in length",
        ));
    }
    let output = rows
        .outer_iter()
        .zip(present)
        .map(|(row, exists)| exists.then_some([row[0], row[1], row[2]]))
        .collect();
    drop(positions);
    Ok(output)
}

fn periodic_error(error: pdbiox::PeriodicError) -> PyErr {
    PyValueError::new_err(format!("invalid periodic data: {error:?}"))
}

fn rotation_error(error: pdbiox::RotationError) -> PyErr {
    PyValueError::new_err(format!("invalid rotation data: {error:?}"))
}

fn eigen_error(error: pdbiox::EigenError) -> PyErr {
    PyValueError::new_err(format!("eigendecomposition failed: {error:?}"))
}

impl From<pdbiox::CircularSummary> for PyCircularSummary {
    fn from(value: pdbiox::CircularSummary) -> Self {
        Self {
            mean: value.mean.map(pdbiox::PeriodicAngle::radians),
            resultant: value.resultant,
            variance: value.variance,
        }
    }
}

impl From<pdbiox::BackboneFrame> for PyBackboneFrame {
    fn from(value: pdbiox::BackboneFrame) -> Self {
        Self {
            tangent: value.tangent,
            normal: value.normal,
            binormal: value.binormal,
            curvature: value.curvature,
            torsion: value.torsion,
        }
    }
}

impl From<pdbiox::HelixGeometry> for PyHelixGeometry {
    fn from(value: pdbiox::HelixGeometry) -> Self {
        Self {
            axis: value.axis,
            rise: value.rise,
            twist: value.twist,
        }
    }
}

impl From<pdbiox::BackboneTorsions> for PyBackboneTorsions {
    fn from(value: pdbiox::BackboneTorsions) -> Self {
        Self {
            phi: value.phi,
            psi: value.psi,
            omega: value.omega,
        }
    }
}
