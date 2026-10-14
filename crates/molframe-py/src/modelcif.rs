//! Predicted aligned error and the `ModelCIF` confidence surface.
//!
//! Everything here is a pass-through: the matrix, its storage choice and every
//! number come from `molframe-modelcif`, and `array()` hands `NumPy` the stored
//! buffer itself rather than a converted copy.

use crate::bindings::PyStructure;
use molframe::formats::modelcif::{
    DensePae, ModelCifExt, PaeKind, PaeMatrix, PaePolicy, PaeStorage, read_pae_json as import_json,
};
use molframe::{AtomAnnotation, MISSING_BAND_CODE, PLDDT_ANNOTATION, PlddtBand};
use numpy::ndarray::ArrayView2;
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyAny;
use std::path::PathBuf;

/// Where the matrix lives, kept alive for as long as any `NumPy` view of it.
#[derive(Debug)]
enum Holder {
    Owned(PaeMatrix),
    Structure(molframe::Structure),
}

/// A square residue-pair aligned-error matrix stored in its narrowest exact layout.
///
/// `array()` is a read-only zero-copy view of the stored cells: `uint8` or
/// `uint16` fixed point (divide by `denominator`), `float16` or `float32`
/// (`denominator` is 1). `decode()` is the only call that allocates a copy.
#[derive(Debug)]
#[pyclass(
    name = "PaeMatrix",
    frozen,
    skip_from_py_object,
    module = "molframe.modelcif"
)]
pub(crate) struct PyPaeMatrix {
    holder: Holder,
}

impl PyPaeMatrix {
    fn matrix(&self) -> PyResult<&PaeMatrix> {
        match &self.holder {
            Holder::Owned(matrix) => Ok(matrix),
            Holder::Structure(structure) => structure
                .model_cif()
                .and_then(|model| model.dense_pae())
                .map(DensePae::matrix)
                .ok_or_else(|| {
                    crate::error::internal("the structure no longer holds a PAE matrix")
                }),
        }
    }
}

fn kind_name(kind: PaeKind) -> &'static str {
    match kind {
        PaeKind::U8 => "uint8",
        PaeKind::U16 => "uint16",
        PaeKind::F16 => "float16",
        PaeKind::F32 => "float32",
    }
}

fn view_of<'py, T: numpy::Element>(
    slf: &Bound<'py, PyPaeMatrix>,
    cells: &[T],
    size: usize,
) -> Bound<'py, PyArray2<T>> {
    // SAFETY: `cells` is `size * size` contiguous elements owned by the holder,
    // which is immutable and kept alive as the array's base object.
    let view = unsafe { ArrayView2::from_shape_ptr((size, size), cells.as_ptr()) };
    // SAFETY: the base object is the holder that owns the allocation.
    let array = unsafe { PyArray2::borrow_from_array(&view, slf.clone().into_any()) };
    array.readwrite().make_nonwriteable();
    array
}

#[pymethods]
impl PyPaeMatrix {
    /// Residues per side.
    #[getter]
    fn size(&self) -> PyResult<usize> {
        Ok(self.matrix()?.size())
    }

    /// `NumPy` dtype name of the stored cells.
    #[getter]
    fn dtype(&self) -> PyResult<&'static str> {
        Ok(kind_name(self.matrix()?.kind()))
    }

    /// Fixed-point divisor of the stored cells; 1 for floating-point storage.
    #[getter]
    fn denominator(&self) -> PyResult<u32> {
        Ok(self.matrix()?.layout().denominator)
    }

    /// Bytes the cell buffer occupies.
    #[getter]
    fn nbytes(&self) -> PyResult<usize> {
        Ok(self.matrix()?.byte_len())
    }

    /// Whether every source value is reproduced exactly.
    #[getter]
    fn is_exact(&self) -> PyResult<bool> {
        Ok(self.matrix()?.is_exact())
    }

    /// Largest absolute difference between a source value and its stored value.
    #[getter]
    fn max_quantisation_error(&self) -> PyResult<f64> {
        Ok(self.matrix()?.max_quantisation_error())
    }

    /// Maximum error the producing model declared, when it stated one.
    #[getter]
    fn declared_maximum(&self) -> PyResult<Option<f64>> {
        Ok(self.matrix()?.declared_maximum())
    }

    /// The stored cells as a read-only `(size, size)` view, with no copy.
    fn array<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyAny>> {
        let this = slf.get();
        let matrix = this.matrix()?;
        let size = matrix.size();
        if let Some(cells) = matrix.as_u8() {
            return Ok(view_of(slf, cells, size).into_any());
        }
        if let Some(cells) = matrix.as_u16() {
            return Ok(view_of(slf, cells, size).into_any());
        }
        if let Some(cells) = matrix.as_f32() {
            return Ok(view_of(slf, cells, size).into_any());
        }
        let bits = matrix
            .as_f16_bits()
            .ok_or_else(|| crate::error::internal("unrecognised PAE storage"))?;
        view_of(slf, bits, size)
            .into_any()
            .call_method1("view", ("float16",))
    }

    /// The value at `(row, column)` as the source stated it.
    fn get(&self, row: usize, column: usize) -> PyResult<f64> {
        self.matrix()?
            .get_f64(row, column)
            .ok_or_else(|| crate::error::index("PAE index out of range"))
    }

    /// Mean of the rectangle `rows x columns` (half-open ranges), or None when empty.
    fn block_mean(&self, rows: (usize, usize), columns: (usize, usize)) -> PyResult<Option<f64>> {
        Ok(self
            .matrix()?
            .block_mean(rows.0..rows.1, columns.0..columns.1))
    }

    /// Mean over every cell.
    fn mean(&self) -> PyResult<Option<f64>> {
        Ok(self.matrix()?.mean())
    }

    /// Decodes row `row` to `float32`; this is a copy.
    fn row(&self, row: usize) -> PyResult<Vec<f32>> {
        let matrix = self.matrix()?;
        let mut out = vec![0.0_f32; matrix.size()];
        matrix
            .row_into(row, &mut out)
            .map_err(|error| crate::error::value(error.to_string()))?;
        Ok(out)
    }

    fn __repr__(&self) -> PyResult<String> {
        let matrix = self.matrix()?;
        Ok(format!(
            "PaeMatrix(size={}, dtype={}, nbytes={})",
            matrix.size(),
            kind_name(matrix.kind()),
            matrix.byte_len()
        ))
    }
}

fn storage(quantise: Option<&str>) -> PyResult<PaeStorage> {
    Ok(match quantise {
        None | Some("exact") => PaeStorage::Exact,
        Some("uint8") => PaeStorage::Quantised(PaeKind::U8),
        Some("uint16") => PaeStorage::Quantised(PaeKind::U16),
        Some("float16") => PaeStorage::Quantised(PaeKind::F16),
        Some("float32") => PaeStorage::Quantised(PaeKind::F32),
        Some(other) => {
            return Err(crate::error::value(format!(
                "quantise must be exact, uint8, uint16, float16 or float32, not {other}"
            )));
        }
    })
}

/// Imports an `AlphaFold`-style PAE JSON file as a dense matrix.
///
/// `quantise` forces a narrower layout (`uint8`, `uint16`, `float16`,
/// `float32`) and accepts the reported error; the default is the narrowest
/// layout that reproduces every value exactly. `max_bytes` bounds the matrix.
#[pyfunction]
#[pyo3(signature = (path, *, quantise=None, max_bytes=None))]
fn read_pae_json(
    py: Python<'_>,
    path: PathBuf,
    quantise: Option<&str>,
    max_bytes: Option<u64>,
) -> PyResult<PyPaeMatrix> {
    let mut policy = PaePolicy::new().with_storage(storage(quantise)?);
    if let Some(limit) = max_bytes {
        policy = policy.with_max_bytes(limit);
    }
    let matrix = py
        .detach(move || {
            let bytes = std::fs::read(&path).map_err(|error| error.to_string())?;
            import_json(&bytes, policy).map_err(|error| error.to_string())
        })
        .map_err(crate::error::value)?;
    Ok(PyPaeMatrix {
        holder: Holder::Owned(matrix),
    })
}

/// The dense PAE matrix a structure read from a `ModelCIF` file carries, or None.
#[pyfunction]
fn pae(structure: &PyStructure) -> Option<PyPaeMatrix> {
    structure
        .inner
        .model_cif()
        .and_then(|model| model.dense_pae())?;
    Some(PyPaeMatrix {
        holder: Holder::Structure(structure.inner.clone()),
    })
}

/// Confidence band of every atom as `uint8`: 0 (`very_low`) to 3 (`very_high`), 255 unscored.
///
/// Returns None when the structure carries no pLDDT. The classification runs in Rust.
#[pyfunction]
fn plddt_bands<'py>(py: Python<'py>, structure: &PyStructure) -> Option<Bound<'py, PyArray1<u8>>> {
    match structure.inner.annotations().get(PLDDT_ANNOTATION) {
        Some(AtomAnnotation::Real(column)) => {
            Some(PyArray1::from_vec(py, PlddtBand::codes(column)))
        }
        _ => None,
    }
}

/// Band names indexed by the codes `plddt_bands` returns.
#[pyfunction]
fn plddt_band_names() -> Vec<&'static str> {
    PlddtBand::ALL.iter().map(|band| band.name()).collect()
}

/// Conventional sRGB colour of each band as a `(4, 3)` `uint8` array.
#[pyfunction]
fn plddt_band_colors(py: Python<'_>) -> PyResult<Bound<'_, PyArray2<u8>>> {
    let flat: Vec<u8> = PlddtBand::ALL.iter().flat_map(|band| band.rgb()).collect();
    PyArray1::from_vec(py, flat).reshape((PlddtBand::ALL.len(), 3))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPaeMatrix>()?;
    module.add_function(wrap_pyfunction!(read_pae_json, module)?)?;
    module.add_function(wrap_pyfunction!(pae, module)?)?;
    module.add_function(wrap_pyfunction!(plddt_bands, module)?)?;
    module.add_function(wrap_pyfunction!(plddt_band_names, module)?)?;
    module.add_function(wrap_pyfunction!(plddt_band_colors, module)?)?;
    module.add("PLDDT_UNSCORED", MISSING_BAND_CODE)
}
