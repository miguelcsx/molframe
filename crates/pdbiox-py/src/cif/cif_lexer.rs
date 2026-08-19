//! Owned, native CIF lexer tokens for Python callers.

use crate::cif_document::PyCifQuoting;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type PositionTuple = (u32, u32, u32);
type SpanTuple = (u32, u32, u32, u32);

#[pyclass(name = "CifToken", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifToken {
    #[pyo3(get)]
    pub(crate) kind: String,
    #[pyo3(get)]
    pub(crate) text: Option<String>,
    #[pyo3(get)]
    pub(crate) quoting: Option<PyCifQuoting>,
}

#[pyclass(name = "Spanned", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifSpanned {
    #[pyo3(get)]
    pub(crate) token: PyCifToken,
    #[pyo3(get)]
    pub(crate) span: SpanTuple,
}

#[derive(Clone, Debug)]
struct OwnedToken {
    spanned: PyCifSpanned,
    position: PositionTuple,
}

#[pyclass(name = "CifLexer")]
pub(crate) struct PyCifLexer {
    tokens: Vec<OwnedToken>,
    cursor: usize,
    position: PositionTuple,
}

#[pymethods]
impl PyCifLexer {
    #[new]
    fn new(data: Vec<u8>) -> PyResult<Self> {
        let mut lexer = pdbiox::cif::lexer::Lexer::new(&data).map_err(lex_error)?;
        let mut tokens = Vec::new();
        loop {
            let Some(spanned) = lexer.next_token().map_err(lex_error)? else {
                break;
            };
            tokens.push(OwnedToken {
                spanned: PyCifSpanned {
                    token: token(spanned.token),
                    span: span(spanned.span),
                },
                position: position(lexer.position()),
            });
        }
        Ok(Self {
            tokens,
            cursor: 0,
            position: (0, 1, 1),
        })
    }

    #[getter]
    fn position(&self) -> PositionTuple {
        self.position
    }

    fn __len__(&self) -> usize {
        self.tokens.len()
    }

    fn next_token(&mut self) -> Option<PyCifSpanned> {
        let token = self.tokens.get(self.cursor)?.clone();
        self.cursor = self.cursor.saturating_add(1);
        self.position = token.position;
        Some(token.spanned)
    }

    fn remaining(&self) -> usize {
        self.tokens.len().saturating_sub(self.cursor)
    }
}

fn token(value: pdbiox::cif::lexer::Token<'_>) -> PyCifToken {
    match value {
        pdbiox::cif::lexer::Token::Block(text) => value_with("block", text, None),
        pdbiox::cif::lexer::Token::FrameStart(text) => value_with("frame_start", text, None),
        pdbiox::cif::lexer::Token::FrameEnd => value_with("frame_end", "", None),
        pdbiox::cif::lexer::Token::Loop => value_with("loop", "", None),
        pdbiox::cif::lexer::Token::Tag(text) => value_with("tag", text, None),
        pdbiox::cif::lexer::Token::Value(text, quoting) => {
            value_with("value", text, Some(quoting.into()))
        }
    }
}

fn value_with(kind: &str, text: &str, quoting: Option<PyCifQuoting>) -> PyCifToken {
    PyCifToken {
        kind: kind.to_owned(),
        text: (!text.is_empty()).then(|| text.to_owned()),
        quoting,
    }
}

fn position(value: pdbiox::Position) -> PositionTuple {
    (value.byte_offset, value.line, value.column)
}

fn span(value: pdbiox::ByteSpan) -> SpanTuple {
    (
        value.start.byte_offset,
        value.start.line,
        value.start.column,
        value.end,
    )
}

fn lex_error(error: pdbiox::cif::lexer::LexError) -> PyErr {
    PyValueError::new_err(format!("{error:?}"))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCifToken>()?;
    module.add_class::<PyCifSpanned>()?;
    module.add_class::<PyCifLexer>()
}
