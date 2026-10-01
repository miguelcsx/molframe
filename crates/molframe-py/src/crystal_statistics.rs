//! Mechanical adapters for resolution shells and amplitude normalization.

use crate::crystal::{PySpaceGroup, PyUnitCell};
use molframe::crystal::{
    BinMethod, ReflectionBinningError, ResolutionBinner, amplitude_normalizers,
};
use numpy::{PyArray1, PyArray2, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
use pyo3::{exceptions::PyValueError, prelude::*};

fn failure(error: ReflectionBinningError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn parse_method(name: &str) -> PyResult<BinMethod> {
    match name {
        "equal_count" => Ok(BinMethod::EqualCount),
        "dstar" => Ok(BinMethod::Dstar),
        "dstar2" => Ok(BinMethod::Dstar2),
        "dstar3" => Ok(BinMethod::Dstar3),
        _ => Err(PyValueError::new_err(
            "method must be one of equal_count, dstar, dstar2, dstar3",
        )),
    }
}

fn rows(hkl: &Bound<'_, PyArray2<i32>>) -> PyResult<Vec<[i32; 3]>> {
    let hkl = hkl.readonly();
    let shape = hkl.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err("hkl must have shape (n, 3)"));
    }
    Ok(hkl
        .as_array()
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

/// Resolution shells over `1/d²` for a set of reflections.
#[derive(Clone, Debug)]
#[pyclass(
    name = "ResolutionBins",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PyResolutionBins(ResolutionBinner);

#[pymethods]
impl PyResolutionBins {
    #[new]
    #[pyo3(signature = (cell, hkl, *, bins=20, method="equal_count"))]
    fn new(
        py: Python<'_>,
        cell: &PyUnitCell,
        hkl: &Bound<'_, PyArray2<i32>>,
        bins: usize,
        method: &str,
    ) -> PyResult<Self> {
        let rows = rows(hkl)?;
        let method = parse_method(method)?;
        let transform = cell.0;
        py.detach(move || ResolutionBinner::for_reflections(method, bins, &transform, &rows))
            .map(Self)
            .map_err(failure)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Upper `1/d²` limit of each shell; the last is infinite.
    #[getter]
    fn limits<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.limits().to_pyarray(py)
    }

    /// Midpoint `1/d²` of each shell.
    #[getter]
    fn midpoints<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.midpoints().to_pyarray(py)
    }

    /// The shell of every `1/d²` value.
    fn indices<'py>(
        &self,
        py: Python<'py>,
        inverse_d2: &Bound<'py, PyArray1<f64>>,
    ) -> PyResult<Bound<'py, PyArray1<i64>>> {
        let readonly = inverse_d2.readonly();
        let slice = readonly.as_slice()?;
        let indices = self
            .0
            .bin_indices(slice)
            .into_iter()
            .map(|bin| {
                i64::try_from(bin).map_err(|_| PyValueError::new_err("shell index overflow"))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(indices.to_pyarray(py))
    }

    /// Smallest `d` in a shell, in ångström.
    fn d_min(&self, bin: usize) -> PyResult<f64> {
        self.0
            .d_min_of_bin(bin)
            .ok_or_else(|| PyValueError::new_err("shell index out of range"))
    }

    /// Largest `d` in a shell, in ångström.
    fn d_max(&self, bin: usize) -> PyResult<f64> {
        self.0
            .d_max_of_bin(bin)
            .ok_or_else(|| PyValueError::new_err("shell index out of range"))
    }
}

/// Multipliers turning amplitudes into normalized amplitudes `E = F · m`.
///
/// Reflections whose amplitude is not finite get `NaN`.
#[pyfunction]
fn normalizers<'py>(
    py: Python<'py>,
    cell: &PyUnitCell,
    space_group: &PySpaceGroup,
    hkl: &Bound<'py, PyArray2<i32>>,
    amplitudes: &Bound<'py, PyArray1<f64>>,
    bins: &PyResolutionBins,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let rows = rows(hkl)?;
    let amplitudes = amplitudes.readonly().as_slice()?.to_vec();
    let (transform, symmetry, binner) = (cell.0, space_group.0.clone(), bins.0.clone());
    let multipliers = py
        .detach(move || amplitude_normalizers(&transform, &symmetry, &rows, &amplitudes, &binner))
        .map_err(failure)?;
    Ok(multipliers.to_pyarray(py))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyResolutionBins>()?;
    module.add_function(wrap_pyfunction!(normalizers, module)?)
}
