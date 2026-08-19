//! Native bindings for path-backed GSD and TNG trajectory containers.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;
use std::path::PathBuf;

pyo3::create_exception!(_native, GsdError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, TngError, pyo3::exceptions::PyValueError);

#[pyclass(name = "GsdOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyGsdOptions {
    #[pyo3(get)]
    pub(crate) length_to_angstrom: f64,
}

#[pymethods]
impl PyGsdOptions {
    #[new]
    fn new(length_to_angstrom: f64) -> PyResult<Self> {
        pdbiox::traj::GsdOptions::new(length_to_angstrom)
            .map(|_| Self { length_to_angstrom })
            .map_err(|error| GsdError::new_err(error.to_string()))
    }
}

impl From<PyGsdOptions> for pdbiox::traj::GsdOptions {
    fn from(value: PyGsdOptions) -> Self {
        Self {
            length_to_angstrom: value.length_to_angstrom,
        }
    }
}

#[pyclass(name = "GsdTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGsdTrajectory {
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
    #[pyo3(get)]
    pub(crate) steps: Vec<u64>,
}

impl TryFrom<pdbiox::traj::GsdTrajectory> for PyGsdTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::GsdTrajectory) -> Result<Self, Self::Error> {
        let frames = value
            .frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            frames,
            steps: value.steps,
        })
    }
}

#[pyfunction]
pub(crate) fn parse_gsd(path: PathBuf, options: PyGsdOptions) -> PyResult<PyGsdTrajectory> {
    pdbiox::traj::parse_gsd(&path, options.into())
        .map_err(|error| GsdError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
pub(crate) fn write_gsd(
    path: PathBuf,
    frames: Vec<PyTimestep>,
    options: PyGsdOptions,
) -> PyResult<()> {
    let frames = frames
        .into_iter()
        .map(TryInto::try_into)
        .collect::<PyResult<Vec<_>>>()?;
    pdbiox::traj::write_gsd(&path, &frames, options.into())
        .map_err(|error| GsdError::new_err(error.to_string()))
}

#[pyclass(name = "TngCompression", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyTngCompression {
    kind: TngCompressionKind,
    precision: f64,
}

#[derive(Clone, Copy, Debug)]
enum TngCompressionKind {
    Uncompressed,
    Lossless,
    Lossy,
}

#[pymethods]
impl PyTngCompression {
    #[staticmethod]
    fn uncompressed() -> Self {
        Self {
            kind: TngCompressionKind::Uncompressed,
            precision: 0.0,
        }
    }

    #[staticmethod]
    fn lossless() -> Self {
        Self {
            kind: TngCompressionKind::Lossless,
            precision: 0.0,
        }
    }

    #[staticmethod]
    fn lossy(precision: f64) -> PyResult<Self> {
        if precision.is_finite() && precision > 0.0 {
            Ok(Self {
                kind: TngCompressionKind::Lossy,
                precision,
            })
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "TNG lossy precision must be finite and positive",
            ))
        }
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.kind {
            TngCompressionKind::Uncompressed => "uncompressed",
            TngCompressionKind::Lossless => "lossless",
            TngCompressionKind::Lossy => "lossy",
        }
    }

    #[getter]
    const fn precision(&self) -> Option<f64> {
        match self.kind {
            TngCompressionKind::Lossy => Some(self.precision),
            TngCompressionKind::Uncompressed | TngCompressionKind::Lossless => None,
        }
    }
}

impl From<PyTngCompression> for pdbiox::traj::TngCompression {
    fn from(value: PyTngCompression) -> Self {
        match value.kind {
            TngCompressionKind::Uncompressed => Self::Uncompressed,
            TngCompressionKind::Lossless => Self::Lossless,
            TngCompressionKind::Lossy => Self::Lossy {
                precision: value.precision,
            },
        }
    }
}

impl From<pdbiox::traj::TngCompression> for PyTngCompression {
    fn from(value: pdbiox::traj::TngCompression) -> Self {
        match value {
            pdbiox::traj::TngCompression::Uncompressed => Self::uncompressed(),
            pdbiox::traj::TngCompression::Lossless => Self::lossless(),
            pdbiox::traj::TngCompression::Lossy { precision } => Self {
                kind: TngCompressionKind::Lossy,
                precision,
            },
        }
    }
}

#[pyclass(name = "TngWriteOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyTngWriteOptions {
    #[pyo3(get)]
    pub(crate) distance_unit_exponent: i64,
    #[pyo3(get)]
    pub(crate) compression: PyTngCompression,
    #[pyo3(get)]
    pub(crate) hashes: bool,
}

#[pymethods]
impl PyTngWriteOptions {
    #[new]
    #[pyo3(signature = (*, distance_unit_exponent=-9, compression=None, hashes=true))]
    fn new(
        distance_unit_exponent: i64,
        compression: Option<PyTngCompression>,
        hashes: bool,
    ) -> Self {
        let defaults = pdbiox::traj::TngWriteOptions::default();
        Self {
            distance_unit_exponent,
            compression: match compression {
                Some(value) => value,
                None => defaults.compression.into(),
            },
            hashes,
        }
    }
}

impl From<PyTngWriteOptions> for pdbiox::traj::TngWriteOptions {
    fn from(value: PyTngWriteOptions) -> Self {
        Self {
            distance_unit_exponent: value.distance_unit_exponent,
            compression: value.compression.into(),
            hashes: value.hashes,
        }
    }
}

#[pyclass(name = "TngTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTngTrajectory {
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
    #[pyo3(get)]
    pub(crate) steps: Vec<i64>,
    #[pyo3(get)]
    pub(crate) distance_unit_exponent: i64,
    #[pyo3(get)]
    pub(crate) compression_precision: f64,
    #[pyo3(get)]
    pub(crate) compression: PyTngCompression,
}

impl TryFrom<pdbiox::traj::TngTrajectory> for PyTngTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::TngTrajectory) -> Result<Self, Self::Error> {
        let frames = value
            .frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            frames,
            steps: value.steps,
            distance_unit_exponent: value.distance_unit_exponent,
            compression_precision: value.compression_precision,
            compression: value.compression.into(),
        })
    }
}

#[pyfunction]
pub(crate) fn parse_tng(path: PathBuf) -> PyResult<PyTngTrajectory> {
    pdbiox::traj::parse_tng(&path)
        .map_err(|error| TngError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
pub(crate) fn write_tng(
    path: PathBuf,
    frames: Vec<PyTimestep>,
    options: PyTngWriteOptions,
) -> PyResult<()> {
    let frames = frames
        .into_iter()
        .map(TryInto::try_into)
        .collect::<PyResult<Vec<_>>>()?;
    pdbiox::traj::write_tng(&path, &frames, options.into())
        .map_err(|error| TngError::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("GsdError", module.py().get_type::<GsdError>())?;
    module.add("TngError", module.py().get_type::<TngError>())?;
    module.add_class::<PyGsdOptions>()?;
    module.add_class::<PyGsdTrajectory>()?;
    module.add_class::<PyTngCompression>()?;
    module.add_class::<PyTngWriteOptions>()?;
    module.add_class::<PyTngTrajectory>()?;
    module.add_function(wrap_pyfunction!(parse_gsd, module)?)?;
    module.add_function(wrap_pyfunction!(write_gsd, module)?)?;
    module.add_function(wrap_pyfunction!(parse_tng, module)?)?;
    module.add_function(wrap_pyfunction!(write_tng, module)?)?;
    Ok(())
}
