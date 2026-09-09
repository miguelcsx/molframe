//! Typed `BinaryCIF` column codecs with zero-copy `NumPy` input.

use numpy::{IntoPyArray, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyList};

type StringArrayParts = (Vec<PyEncoding>, String, Vec<PyEncoding>, Vec<u8>);

#[pyclass(name = "DataType", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyDataType {
    Int8,
    Int16,
    Int32,
    Uint8,
    Uint16,
    Uint32,
    Float32,
    Float64,
}

impl From<PyDataType> for pdbiox::bcif::DataType {
    fn from(value: PyDataType) -> Self {
        match value {
            PyDataType::Int8 => Self::Int8,
            PyDataType::Int16 => Self::Int16,
            PyDataType::Int32 => Self::Int32,
            PyDataType::Uint8 => Self::Uint8,
            PyDataType::Uint16 => Self::Uint16,
            PyDataType::Uint32 => Self::Uint32,
            PyDataType::Float32 => Self::Float32,
            PyDataType::Float64 => Self::Float64,
        }
    }
}

impl From<pdbiox::bcif::DataType> for PyDataType {
    fn from(value: pdbiox::bcif::DataType) -> Self {
        match value {
            pdbiox::bcif::DataType::Int8 => Self::Int8,
            pdbiox::bcif::DataType::Int16 => Self::Int16,
            pdbiox::bcif::DataType::Int32 => Self::Int32,
            pdbiox::bcif::DataType::Uint8 => Self::Uint8,
            pdbiox::bcif::DataType::Uint16 => Self::Uint16,
            pdbiox::bcif::DataType::Uint32 => Self::Uint32,
            pdbiox::bcif::DataType::Float32 => Self::Float32,
            pdbiox::bcif::DataType::Float64 => Self::Float64,
        }
    }
}

#[pyclass(name = "Encoding", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEncoding(pub(crate) pdbiox::bcif::Encoding);

#[pymethods]
impl PyEncoding {
    #[staticmethod]
    fn byte_array(data_type: PyDataType) -> Self {
        Self(pdbiox::bcif::Encoding::ByteArray {
            r#type: data_type.into(),
        })
    }

    #[staticmethod]
    fn fixed_point(factor: f64, src_type: PyDataType) -> Self {
        Self(pdbiox::bcif::Encoding::FixedPoint {
            factor,
            src_type: src_type.into(),
        })
    }

    #[staticmethod]
    fn interval_quantization(
        minimum: f64,
        maximum: f64,
        num_steps: u32,
        src_type: PyDataType,
    ) -> Self {
        Self(pdbiox::bcif::Encoding::IntervalQuantization {
            min: minimum,
            max: maximum,
            num_steps,
            src_type: src_type.into(),
        })
    }

    #[staticmethod]
    fn run_length(src_type: PyDataType, src_size: usize) -> Self {
        Self(pdbiox::bcif::Encoding::RunLength {
            src_type: src_type.into(),
            src_size,
        })
    }

    #[staticmethod]
    fn delta(origin: i64, src_type: PyDataType) -> Self {
        Self(pdbiox::bcif::Encoding::Delta {
            origin,
            src_type: src_type.into(),
        })
    }

    #[staticmethod]
    fn integer_packing(byte_count: u8, is_unsigned: bool, src_size: usize) -> Self {
        Self(pdbiox::bcif::Encoding::IntegerPacking {
            byte_count,
            is_unsigned,
            src_size,
        })
    }

    #[staticmethod]
    fn string_array(
        data_encoding: Vec<PyEncoding>,
        string_data: String,
        offset_encoding: Vec<PyEncoding>,
        offsets: &Bound<'_, PyBytes>,
    ) -> Self {
        Self(pdbiox::bcif::Encoding::StringArray {
            data_encoding: data_encoding.into_iter().map(|value| value.0).collect(),
            string_data,
            offset_encoding: offset_encoding.into_iter().map(|value| value.0).collect(),
            offsets: offsets.as_bytes().to_vec(),
        })
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.0 {
            pdbiox::bcif::Encoding::ByteArray { .. } => "byte_array",
            pdbiox::bcif::Encoding::FixedPoint { .. } => "fixed_point",
            pdbiox::bcif::Encoding::IntervalQuantization { .. } => "interval_quantization",
            pdbiox::bcif::Encoding::RunLength { .. } => "run_length",
            pdbiox::bcif::Encoding::Delta { .. } => "delta",
            pdbiox::bcif::Encoding::IntegerPacking { .. } => "integer_packing",
            pdbiox::bcif::Encoding::StringArray { .. } => "string_array",
        }
    }

    #[getter]
    fn data_type(&self) -> Option<PyDataType> {
        match self.0 {
            pdbiox::bcif::Encoding::ByteArray { r#type } => Some(r#type.into()),
            _ => None,
        }
    }

    #[getter]
    fn src_type(&self) -> Option<PyDataType> {
        match self.0 {
            pdbiox::bcif::Encoding::FixedPoint { src_type, .. }
            | pdbiox::bcif::Encoding::IntervalQuantization { src_type, .. }
            | pdbiox::bcif::Encoding::RunLength { src_type, .. }
            | pdbiox::bcif::Encoding::Delta { src_type, .. } => Some(src_type.into()),
            _ => None,
        }
    }

    #[getter]
    fn factor(&self) -> Option<f64> {
        match self.0 {
            pdbiox::bcif::Encoding::FixedPoint { factor, .. } => Some(factor),
            _ => None,
        }
    }

    #[getter]
    fn interval(&self) -> Option<(f64, f64, u32)> {
        match self.0 {
            pdbiox::bcif::Encoding::IntervalQuantization {
                min,
                max,
                num_steps,
                ..
            } => Some((min, max, num_steps)),
            _ => None,
        }
    }

    #[getter]
    fn origin(&self) -> Option<i64> {
        match self.0 {
            pdbiox::bcif::Encoding::Delta { origin, .. } => Some(origin),
            _ => None,
        }
    }

    #[getter]
    fn src_size(&self) -> Option<usize> {
        match self.0 {
            pdbiox::bcif::Encoding::RunLength { src_size, .. }
            | pdbiox::bcif::Encoding::IntegerPacking { src_size, .. } => Some(src_size),
            _ => None,
        }
    }

    #[getter]
    fn packing(&self) -> Option<(u8, bool)> {
        match self.0 {
            pdbiox::bcif::Encoding::IntegerPacking {
                byte_count,
                is_unsigned,
                ..
            } => Some((byte_count, is_unsigned)),
            _ => None,
        }
    }

    fn string_array_parts(&self) -> Option<StringArrayParts> {
        match &self.0 {
            pdbiox::bcif::Encoding::StringArray {
                data_encoding,
                string_data,
                offset_encoding,
                offsets,
            } => Some((
                data_encoding.iter().cloned().map(PyEncoding).collect(),
                string_data.clone(),
                offset_encoding.iter().cloned().map(PyEncoding).collect(),
                offsets.clone(),
            )),
            _ => None,
        }
    }
}

#[pyclass(name = "EncodedData", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEncodedData(pub(crate) pdbiox::bcif::EncodedData);

#[pymethods]
impl PyEncodedData {
    #[new]
    fn new(encoding: Vec<PyEncoding>, data: &Bound<'_, PyBytes>) -> Self {
        Self(pdbiox::bcif::EncodedData {
            encoding: encoding.into_iter().map(|value| value.0).collect(),
            data: data.as_bytes().to_vec(),
        })
    }

    #[getter]
    fn encoding(&self) -> Vec<PyEncoding> {
        self.0.encoding.iter().cloned().map(PyEncoding).collect()
    }

    #[getter]
    fn data<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.0.data)
    }
}

#[pyclass(name = "Decoded", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyDecoded {
    #[pyo3(get)]
    kind: &'static str,
    pub(crate) values: Py<PyAny>,
}

#[pymethods]
impl PyDecoded {
    #[getter]
    fn values(&self, py: Python<'_>) -> Py<PyAny> {
        self.values.clone_ref(py)
    }
}

#[pyfunction]
pub(crate) fn encode_integers(
    py: Python<'_>,
    values: PyReadonlyArray1<'_, i64>,
) -> PyResult<PyEncodedData> {
    let values = contiguous(&values)?;
    py.detach(|| pdbiox::bcif::encode_integers(values))
        .map(PyEncodedData)
        .map_err(codec_error)
}

#[pyfunction]
pub(crate) fn encode_floats(
    py: Python<'_>,
    values: PyReadonlyArray1<'_, f64>,
) -> PyResult<PyEncodedData> {
    let values = contiguous(&values)?;
    py.detach(|| pdbiox::bcif::encode_floats(values))
        .map(PyEncodedData)
        .map_err(codec_error)
}

#[pyfunction]
pub(crate) fn encode_interval(
    py: Python<'_>,
    values: PyReadonlyArray1<'_, f64>,
    steps: u32,
) -> PyResult<PyEncodedData> {
    let values = contiguous(&values)?;
    py.detach(|| pdbiox::bcif::encode_interval(values, steps))
        .map(PyEncodedData)
        .map_err(codec_error)
}

#[pyfunction]
pub(crate) fn encode_strings(py: Python<'_>, values: Vec<String>) -> PyResult<PyEncodedData> {
    py.detach(move || pdbiox::bcif::encode_strings(&values))
        .map(PyEncodedData)
        .map_err(codec_error)
}

#[pyfunction]
pub(crate) fn decode(py: Python<'_>, encoded: &PyEncodedData) -> PyResult<PyDecoded> {
    let encoded = encoded.0.clone();
    py.detach(move || pdbiox::bcif::decode(&encoded))
        .map_err(codec_error)
        .and_then(|decoded| decoded_into_python(py, decoded))
}

fn contiguous<'a, T>(array: &'a PyReadonlyArray1<'_, T>) -> PyResult<&'a [T]>
where
    T: numpy::Element,
{
    array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "values must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })
}

fn decoded_into_python(py: Python<'_>, decoded: pdbiox::bcif::Decoded) -> PyResult<PyDecoded> {
    let (kind, values) = match decoded {
        pdbiox::bcif::Decoded::Integers(values) => {
            ("integers", values.into_pyarray(py).into_any().unbind())
        }
        pdbiox::bcif::Decoded::Floats(values) => {
            ("floats", values.into_pyarray(py).into_any().unbind())
        }
        pdbiox::bcif::Decoded::Strings(values) => ("strings", strings_into_python(py, &values)?),
    };
    Ok(PyDecoded { kind, values })
}

fn strings_into_python(
    py: Python<'_>,
    values: &pdbiox::bcif::DecodedStringColumn,
) -> PyResult<Py<PyAny>> {
    let mut rows = Vec::with_capacity(values.indices().len());
    for index in values.indices() {
        let Ok(index) = usize::try_from(*index) else {
            return Err(PyValueError::new_err(
                "decoded string index exceeds the platform address space",
            ));
        };
        let Some(value) = values.dictionary().get(index) else {
            return Err(PyValueError::new_err(
                "decoded string index is outside the dictionary",
            ));
        };
        rows.push(value.as_ref());
    }
    Ok(PyList::new(py, rows)?.into_any().unbind())
}

fn codec_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDataType>()?;
    module.add_class::<PyEncoding>()?;
    module.add_class::<PyEncodedData>()?;
    module.add_class::<PyDecoded>()?;
    module.add_function(wrap_pyfunction!(decode, module)?)?;
    module.add_function(wrap_pyfunction!(encode_integers, module)?)?;
    module.add_function(wrap_pyfunction!(encode_floats, module)?)?;
    module.add_function(wrap_pyfunction!(encode_interval, module)?)?;
    module.add_function(wrap_pyfunction!(encode_strings, module)?)?;
    Ok(())
}
