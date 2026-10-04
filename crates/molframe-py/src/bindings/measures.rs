//! Measurements and fits over `NumPy` coordinate arrays.
//!
//! Every function takes and returns arrays; angles are in degrees, and a
//! quantity that does not exist for its input (an angle at a coincident
//! vertex, a torsion through a collinear triple) is `NaN`, never a placeholder.

use super::geometry::coordinates;
use molframe::geometry as geom;
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyUntypedArrayMethods,
    ToPyArray,
};
use pyo3::prelude::*;

type Coordinates<'py> = Bound<'py, PyArray2<f32>>;

fn masses_or_ones(masses: Option<PyReadonlyArray1<'_, f64>>, count: usize) -> PyResult<Vec<f64>> {
    match masses {
        Some(values) => {
            let values = values.as_slice().map_err(|_| {
                crate::error::value("masses must be C-contiguous; call numpy.ascontiguousarray")
            })?;
            if values.len() == count {
                Ok(values.to_vec())
            } else {
                Err(crate::error::value(
                    "masses must have one entry per position",
                ))
            }
        }
        None => Ok(vec![1.0; count]),
    }
}

fn to_degrees(radians: Vec<f64>) -> Vec<f64> {
    radians.into_iter().map(geom::degrees).collect()
}

/// Distances between corresponding rows of two `(n, 3)` arrays, in ångström.
#[pyfunction]
pub(crate) fn distances<'py>(
    py: Python<'py>,
    left: &Coordinates<'py>,
    right: &Coordinates<'py>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let (left, right) = (left.readonly(), right.readonly());
    let (left, right) = (coordinates(&left)?, coordinates(&right)?);
    let mut output = vec![0.0; left.len()];
    py.detach(|| geom::distances_into(left, right, &mut output))
        .map_err(crate::error::kernel)?;
    Ok(output.into_pyarray(py))
}

/// Angles at `vertex` between rays to `first` and `third`, in degrees.
#[pyfunction]
pub(crate) fn angles<'py>(
    py: Python<'py>,
    first: &Coordinates<'py>,
    vertex: &Coordinates<'py>,
    third: &Coordinates<'py>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let (first, vertex, third) = (first.readonly(), vertex.readonly(), third.readonly());
    let (first, vertex, third) = (
        coordinates(&first)?,
        coordinates(&vertex)?,
        coordinates(&third)?,
    );
    let mut output = vec![0.0; first.len()];
    py.detach(|| geom::angles_into(first, vertex, third, &mut output))
        .map_err(crate::error::kernel)?;
    Ok(to_degrees(output).into_pyarray(py))
}

/// Signed torsions through four rows of points, in degrees, in `(-180, 180]`.
#[pyfunction]
pub(crate) fn dihedrals<'py>(
    py: Python<'py>,
    first: &Coordinates<'py>,
    second: &Coordinates<'py>,
    third: &Coordinates<'py>,
    fourth: &Coordinates<'py>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let (a, b, c, d) = (
        first.readonly(),
        second.readonly(),
        third.readonly(),
        fourth.readonly(),
    );
    let (a, b, c, d) = (
        coordinates(&a)?,
        coordinates(&b)?,
        coordinates(&c)?,
        coordinates(&d)?,
    );
    let mut output = vec![0.0; a.len()];
    py.detach(|| geom::torsions_into(a, b, c, d, &mut output))
        .map_err(crate::error::kernel)?;
    Ok(to_degrees(output).into_pyarray(py))
}

/// The mass-weighted centre; `masses=None` weights every position equally.
#[pyfunction]
#[pyo3(signature = (positions, masses=None))]
pub(crate) fn centre_of_mass(
    py: Python<'_>,
    positions: &Coordinates<'_>,
    masses: Option<PyReadonlyArray1<'_, f64>>,
) -> PyResult<Option<[f64; 3]>> {
    let positions = positions.readonly();
    let positions = coordinates(&positions)?;
    let masses = masses_or_ones(masses, positions.len())?;
    Ok(py.detach(|| geom::centre_of_mass(positions, &masses)))
}

/// The mass-weighted radius of gyration in ångström.
#[pyfunction]
#[pyo3(signature = (positions, masses=None))]
pub(crate) fn radius_of_gyration(
    py: Python<'_>,
    positions: &Coordinates<'_>,
    masses: Option<PyReadonlyArray1<'_, f64>>,
) -> PyResult<Option<f64>> {
    let positions = positions.readonly();
    let positions = coordinates(&positions)?;
    let masses = masses_or_ones(masses, positions.len())?;
    Ok(py.detach(|| geom::radius_of_gyration(positions, &masses)))
}

/// The `(3, 3)` inertia tensor about the centre of mass.
#[pyfunction]
#[pyo3(signature = (positions, masses=None))]
pub(crate) fn inertia_tensor<'py>(
    py: Python<'py>,
    positions: &Coordinates<'py>,
    masses: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<Option<Bound<'py, PyArray2<f64>>>> {
    let coordinates_view = positions.readonly();
    let positions = coordinates(&coordinates_view)?;
    let masses = masses_or_ones(masses, positions.len())?;
    let tensor = py.detach(|| geom::inertia_tensor(positions, &masses));
    tensor
        .map(|rows| {
            let flat: Vec<f64> = rows.into_iter().flatten().collect();
            PyArray1::from_vec(py, flat).reshape([3, 3])
        })
        .transpose()
}

/// The principal moments and axes: `(moments, axes)` with the axes as columns,
/// longest first.
#[pyfunction]
#[pyo3(signature = (positions, masses=None))]
pub(crate) fn principal_axes<'py>(
    py: Python<'py>,
    positions: &Coordinates<'py>,
    masses: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<Axes<'py>> {
    let view = positions.readonly();
    let positions = coordinates(&view)?;
    let masses = masses_or_ones(masses, positions.len())?;
    let found = py
        .detach(|| geom::principal_axes(positions, &masses))
        .map_err(crate::error::kernel)?;
    decomposition(py, found)
}

/// The gyration tensor's eigenvalues and axes: `(values, axes)`.
#[pyfunction]
pub(crate) fn gyration_axes<'py>(
    py: Python<'py>,
    positions: &Coordinates<'py>,
) -> PyResult<Axes<'py>> {
    let view = positions.readonly();
    let positions = coordinates(&view)?;
    let found = py
        .detach(|| geom::gyration_axes(positions))
        .map_err(crate::error::kernel)?;
    decomposition(py, found)
}

type Axes<'py> = Option<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray2<f64>>)>;

fn decomposition(py: Python<'_>, found: Option<geom::Decomposition<3>>) -> PyResult<Axes<'_>> {
    let Some(found) = found else {
        return Ok(None);
    };
    let flat: Vec<f64> = found.vectors.into_iter().flatten().collect();
    Ok(Some((
        found.values.to_pyarray(py),
        PyArray1::from_vec(py, flat).reshape([3, 3])?,
    )))
}

/// How far the shape departs from a sphere, `0` (spherical) to `1` (a line).
#[pyfunction]
pub(crate) fn asphericity(py: Python<'_>, positions: &Coordinates<'_>) -> PyResult<Option<f64>> {
    let view = positions.readonly();
    let positions = coordinates(&view)?;
    py.detach(|| geom::asphericity(positions))
        .map_err(crate::error::kernel)
}

/// The shape parameter of the gyration tensor: negative for oblate shapes,
/// positive for prolate ones.
#[pyfunction]
pub(crate) fn shape_parameter(
    py: Python<'_>,
    positions: &Coordinates<'_>,
) -> PyResult<Option<f64>> {
    let view = positions.readonly();
    let positions = coordinates(&view)?;
    py.detach(|| geom::shape_parameter(positions))
        .map_err(crate::error::kernel)
}

/// The best-fit plane through at least three points: `(centre, normal)`.
#[pyfunction]
pub(crate) fn best_fit_plane(
    py: Python<'_>,
    points: &Coordinates<'_>,
) -> PyResult<Option<([f64; 3], [f64; 3])>> {
    let view = points.readonly();
    let points = coordinates(&view)?;
    let plane = py
        .detach(|| geom::best_fit_plane(points))
        .map_err(crate::error::kernel)?;
    Ok(plane.map(|plane| (plane.centre, plane.normal)))
}

/// The RMS distance of the points from their best-fit plane, in ångström.
#[pyfunction]
pub(crate) fn plane_deviation(py: Python<'_>, points: &Coordinates<'_>) -> PyResult<Option<f64>> {
    let view = points.readonly();
    let points = coordinates(&view)?;
    py.detach(|| geom::plane_deviation(points))
        .map_err(crate::error::kernel)
}

/// The root-mean-square deviation after the optimal rigid superposition.
#[pyfunction]
pub(crate) fn rmsd_after_fit(
    py: Python<'_>,
    mobile: &Coordinates<'_>,
    reference: &Coordinates<'_>,
) -> PyResult<f64> {
    let (mobile, reference) = (mobile.readonly(), reference.readonly());
    let (mobile, reference) = (coordinates(&mobile)?, coordinates(&reference)?);
    py.detach(|| geom::rmsd_after_fit(mobile, reference))
        .map_err(crate::error::kernel)
}

/// A rigid transform and the deviation it achieves.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Superposition",
    frozen,
    skip_from_py_object,
    module = "molframe.geometry"
)]
pub(crate) struct PySuperposition {
    fit: geom::Superposition,
}

#[pymethods]
impl PySuperposition {
    /// The `(3, 3)` rotation, applied to a column vector: `R @ x + t`.
    #[getter]
    fn rotation<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let flat: Vec<f64> = self.fit.transform.rotation.into_iter().flatten().collect();
        PyArray1::from_vec(py, flat).reshape([3, 3])
    }

    /// The translation applied after rotating.
    #[getter]
    const fn translation(&self) -> [f64; 3] {
        self.fit.transform.translation
    }

    /// The deviation after the transform is applied, in ångström.
    #[getter]
    const fn rmsd(&self) -> f64 {
        self.fit.rmsd
    }

    /// The transformed copy of an `(n, 3)` array.
    fn apply<'py>(
        &self,
        py: Python<'py>,
        positions: &Coordinates<'py>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let view = positions.readonly();
        let mut moved = coordinates(&view)?.to_vec();
        let transform = self.fit.transform;
        py.detach(|| transform.apply_all(&mut moved));
        let rows = moved.len();
        let flat: Vec<f32> = moved.into_iter().flatten().collect();
        PyArray1::from_vec(py, flat).reshape([rows, 3])
    }

    fn __repr__(&self) -> String {
        format!("Superposition(rmsd={})", self.fit.rmsd)
    }
}

/// The rigid transform that best carries `mobile` onto `reference`.
///
/// The two sets must already correspond row by row; which atom matches which
/// is a separate question that geometry cannot answer.
#[pyfunction]
pub(crate) fn superpose(
    py: Python<'_>,
    mobile: &Coordinates<'_>,
    reference: &Coordinates<'_>,
) -> PyResult<PySuperposition> {
    let (mobile, reference) = (mobile.readonly(), reference.readonly());
    let (mobile, reference) = (coordinates(&mobile)?, coordinates(&reference)?);
    py.detach(|| geom::superpose(mobile, reference))
        .map(|fit| PySuperposition { fit })
        .map_err(crate::error::kernel)
}

/// Per-atom root-mean-square fluctuation about the mean over `(frames, atoms, 3)`
/// positions, in ångström. Frames are used as given; superpose them first to
/// remove rigid motion.
#[pyfunction]
pub(crate) fn rmsf<'py>(
    py: Python<'py>,
    positions: &Bound<'py, numpy::PyArray3<f32>>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let view = positions.readonly();
    let shape = view.shape().to_vec();
    let flat = view.as_slice().map_err(|_| {
        crate::error::value("positions must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    let atoms = shape[1];
    let frames: Vec<&[[f32; 3]]> = if atoms == 0 {
        vec![&[]; shape[0]]
    } else {
        flat.as_chunks::<3>().0.chunks(atoms).collect()
    };
    let values = py
        .detach(|| geom::rmsf(&frames))
        .map_err(crate::error::kernel)?;
    Ok(values.into_pyarray(py))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySuperposition>()?;
    module.add_function(wrap_pyfunction!(distances, module)?)?;
    module.add_function(wrap_pyfunction!(angles, module)?)?;
    module.add_function(wrap_pyfunction!(dihedrals, module)?)?;
    module.add_function(wrap_pyfunction!(centre_of_mass, module)?)?;
    module.add_function(wrap_pyfunction!(radius_of_gyration, module)?)?;
    module.add_function(wrap_pyfunction!(inertia_tensor, module)?)?;
    module.add_function(wrap_pyfunction!(principal_axes, module)?)?;
    module.add_function(wrap_pyfunction!(gyration_axes, module)?)?;
    module.add_function(wrap_pyfunction!(asphericity, module)?)?;
    module.add_function(wrap_pyfunction!(shape_parameter, module)?)?;
    module.add_function(wrap_pyfunction!(best_fit_plane, module)?)?;
    module.add_function(wrap_pyfunction!(plane_deviation, module)?)?;
    module.add_function(wrap_pyfunction!(rmsd_after_fit, module)?)?;
    module.add_function(wrap_pyfunction!(superpose, module)?)?;
    module.add_function(wrap_pyfunction!(rmsf, module)?)
}
