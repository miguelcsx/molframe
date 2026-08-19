//! Typed alphabet and validated sequence adapters.

use ::pdbiox;
use ::pdbiox::seq::Alphabet;
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBytes, PyInt};

create_exception!(_native, AlphabetError, PyValueError);
create_exception!(_native, SequenceError, PyValueError);

fn bytes<'py>(py: Python<'py>, value: &[u8]) -> Bound<'py, PyBytes> {
    PyBytes::new(py, value)
}

macro_rules! static_alphabet {
    ($rust:ident, $python:ident, $name:literal) => {
        #[pyclass(name = $name, frozen, from_py_object)]
        #[derive(Clone, Copy, Debug)]
        pub(crate) struct $python(pub(crate) pdbiox::seq::$rust);

        #[pymethods]
        impl $python {
            #[new]
            fn new() -> Self {
                Self(pdbiox::seq::$rust)
            }

            fn symbols<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
                bytes(
                    py,
                    <pdbiox::seq::$rust as pdbiox::seq::Alphabet>::symbols(&self.0),
                )
            }

            fn encode_symbol(&self, symbol: u8) -> Option<u8> {
                <pdbiox::seq::$rust as pdbiox::seq::Alphabet>::encode_symbol(&self.0, symbol)
            }

            fn decode_symbol(&self, code: u8) -> Option<u8> {
                <pdbiox::seq::$rust as pdbiox::seq::Alphabet>::decode_symbol(&self.0, code)
            }

            fn encode<'py>(
                &self,
                py: Python<'py>,
                symbols: Vec<u8>,
            ) -> Option<Bound<'py, PyBytes>> {
                <pdbiox::seq::$rust as pdbiox::seq::Alphabet>::encode(&self.0, &symbols)
                    .as_deref()
                    .map(|value| bytes(py, value))
            }

            fn decode<'py>(&self, py: Python<'py>, codes: Vec<u8>) -> Option<Bound<'py, PyBytes>> {
                <pdbiox::seq::$rust as pdbiox::seq::Alphabet>::decode(&self.0, &codes)
                    .as_deref()
                    .map(|value| bytes(py, value))
            }

            fn len(&self) -> usize {
                <pdbiox::seq::$rust as pdbiox::seq::Alphabet>::len(&self.0)
            }

            fn is_empty(&self) -> bool {
                <pdbiox::seq::$rust as pdbiox::seq::Alphabet>::is_empty(&self.0)
            }

            fn __repr__(&self) -> &'static str {
                concat!($name, "()")
            }
        }
    };
}

static_alphabet!(ProteinAlphabet, PyProteinAlphabet, "ProteinAlphabet");
static_alphabet!(DnaAlphabet, PyDnaAlphabet, "DnaAlphabet");
static_alphabet!(RnaAlphabet, PyRnaAlphabet, "RnaAlphabet");
static_alphabet!(
    NucleotideAlphabet,
    PyNucleotideAlphabet,
    "NucleotideAlphabet"
);

#[pyclass(name = "CustomAlphabet", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCustomAlphabet(pub(crate) pdbiox::seq::CustomAlphabet);

#[pymethods]
impl PyCustomAlphabet {
    #[new]
    fn new(symbols: Vec<u8>) -> PyResult<Self> {
        pdbiox::seq::CustomAlphabet::new(symbols.into_boxed_slice())
            .map(Self)
            .map_err(|error| AlphabetError::new_err(error.to_string()))
    }

    fn symbols<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        bytes(
            py,
            <pdbiox::seq::CustomAlphabet as pdbiox::seq::Alphabet>::symbols(&self.0),
        )
    }

    fn encode_symbol(&self, symbol: u8) -> Option<u8> {
        pdbiox::seq::Alphabet::encode_symbol(&self.0, symbol)
    }

    fn decode_symbol(&self, code: u8) -> Option<u8> {
        pdbiox::seq::Alphabet::decode_symbol(&self.0, code)
    }

    fn encode<'py>(&self, py: Python<'py>, symbols: Vec<u8>) -> Option<Bound<'py, PyBytes>> {
        pdbiox::seq::Alphabet::encode(&self.0, &symbols)
            .as_deref()
            .map(|value| bytes(py, value))
    }

    fn decode<'py>(&self, py: Python<'py>, codes: Vec<u8>) -> Option<Bound<'py, PyBytes>> {
        pdbiox::seq::Alphabet::decode(&self.0, &codes)
            .as_deref()
            .map(|value| bytes(py, value))
    }

    fn len(&self) -> usize {
        pdbiox::seq::Alphabet::len(&self.0)
    }

    fn is_empty(&self) -> bool {
        pdbiox::seq::Alphabet::is_empty(&self.0)
    }
}

#[derive(Clone, Debug)]
enum AlphabetValue {
    Protein(pdbiox::seq::ProteinAlphabet),
    Dna(pdbiox::seq::DnaAlphabet),
    Rna(pdbiox::seq::RnaAlphabet),
    Nucleotide(pdbiox::seq::NucleotideAlphabet),
    Custom(pdbiox::seq::CustomAlphabet),
}

impl AlphabetValue {
    fn encode(&self, symbols: &[u8]) -> Option<Vec<u8>> {
        match self {
            Self::Protein(value) => value.encode(symbols),
            Self::Dna(value) => value.encode(symbols),
            Self::Rna(value) => value.encode(symbols),
            Self::Nucleotide(value) => value.encode(symbols),
            Self::Custom(value) => value.encode(symbols),
        }
    }

    fn decode(&self, codes: &[u8]) -> Option<Vec<u8>> {
        match self {
            Self::Protein(value) => value.decode(codes),
            Self::Dna(value) => value.decode(codes),
            Self::Rna(value) => value.decode(codes),
            Self::Nucleotide(value) => value.decode(codes),
            Self::Custom(value) => value.decode(codes),
        }
    }

    fn name(&self) -> &'static str {
        match self {
            Self::Protein(_) => "ProteinAlphabet",
            Self::Dna(_) => "DnaAlphabet",
            Self::Rna(_) => "RnaAlphabet",
            Self::Nucleotide(_) => "NucleotideAlphabet",
            Self::Custom(_) => "CustomAlphabet",
        }
    }
}

fn extract_alphabet(value: &Bound<'_, PyAny>) -> PyResult<AlphabetValue> {
    if let Ok(value) = value.extract::<PyRef<'_, PyProteinAlphabet>>() {
        return Ok(AlphabetValue::Protein(value.0));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyDnaAlphabet>>() {
        return Ok(AlphabetValue::Dna(value.0));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyRnaAlphabet>>() {
        return Ok(AlphabetValue::Rna(value.0));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyNucleotideAlphabet>>() {
        return Ok(AlphabetValue::Nucleotide(value.0));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyCustomAlphabet>>() {
        return Ok(AlphabetValue::Custom(value.0.clone()));
    }
    Err(pyo3::exceptions::PyTypeError::new_err(
        "alphabet must be a built-in or CustomAlphabet",
    ))
}

#[pyclass(name = "Sequence", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySequence {
    alphabet: AlphabetValue,
    codes: Vec<u8>,
}

#[pymethods]
impl PySequence {
    #[new]
    fn new(alphabet: &Bound<'_, PyAny>, symbols: Vec<u8>) -> PyResult<Self> {
        let alphabet = extract_alphabet(alphabet)?;
        let codes = alphabet.encode(&symbols).ok_or_else(|| {
            SequenceError::new_err("sequence contains a symbol outside its alphabet")
        })?;
        Ok(Self { alphabet, codes })
    }

    #[staticmethod]
    fn from_codes(alphabet: &Bound<'_, PyAny>, codes: Vec<u8>) -> PyResult<Self> {
        let alphabet = extract_alphabet(alphabet)?;
        if alphabet.decode(&codes).is_none() {
            return Err(SequenceError::new_err(
                "sequence contains a code outside its alphabet",
            ));
        }
        Ok(Self { alphabet, codes })
    }

    #[getter]
    fn codes<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        bytes(py, &self.codes)
    }

    #[getter]
    fn symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        self.alphabet
            .decode(&self.codes)
            .map(|value| bytes(py, &value))
            .ok_or_else(|| SequenceError::new_err("sequence storage is invalid"))
    }

    #[getter]
    fn alphabet<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        match &self.alphabet {
            AlphabetValue::Protein(value) => Ok(Py::new(py, PyProteinAlphabet(*value))?
                .into_bound(py)
                .into_any()),
            AlphabetValue::Dna(value) => Ok(Py::new(py, PyDnaAlphabet(*value))?
                .into_bound(py)
                .into_any()),
            AlphabetValue::Rna(value) => Ok(Py::new(py, PyRnaAlphabet(*value))?
                .into_bound(py)
                .into_any()),
            AlphabetValue::Nucleotide(value) => Ok(Py::new(py, PyNucleotideAlphabet(*value))?
                .into_bound(py)
                .into_any()),
            AlphabetValue::Custom(value) => Ok(Py::new(py, PyCustomAlphabet(value.clone()))?
                .into_bound(py)
                .into_any()),
        }
    }

    fn len(&self) -> usize {
        self.codes.len()
    }

    fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }

    fn __repr__(&self) -> String {
        format!(
            "Sequence(alphabet={}, length={})",
            self.alphabet.name(),
            self.codes.len()
        )
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("AlphabetError", module.py().get_type::<AlphabetError>())?;
    module.add("SequenceError", module.py().get_type::<SequenceError>())?;
    module.add_class::<PyProteinAlphabet>()?;
    module.add_class::<PyDnaAlphabet>()?;
    module.add_class::<PyRnaAlphabet>()?;
    module.add_class::<PyNucleotideAlphabet>()?;
    module.add_class::<PyCustomAlphabet>()?;
    module.add_class::<PySequence>()?;
    module.add("Alphabet", module.py().get_type::<PyProteinAlphabet>())?;
    module.add("SequenceCode", module.py().get_type::<PyInt>())?;
    module.add(
        "PROTEIN",
        Py::new(module.py(), PyProteinAlphabet(pdbiox::seq::ProteinAlphabet))?,
    )?;
    module.add(
        "DNA",
        Py::new(module.py(), PyDnaAlphabet(pdbiox::seq::DnaAlphabet))?,
    )?;
    module.add(
        "RNA",
        Py::new(module.py(), PyRnaAlphabet(pdbiox::seq::RnaAlphabet))?,
    )?;
    module.add(
        "NUCLEOTIDE",
        Py::new(
            module.py(),
            PyNucleotideAlphabet(pdbiox::seq::NucleotideAlphabet),
        )?,
    )?;
    Ok(())
}
