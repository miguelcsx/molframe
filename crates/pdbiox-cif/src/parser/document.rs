//! Assembling tokens into a document.
//!
//! Two shapes carry data: a tag followed by one value, and a loop header
//! followed by tags and then rows of values. Everything else is structure.
//!
//! A row count that is not a whole multiple of the item count is the single most
//! common way a hand-edited file breaks, almost always because an unquoted value
//! contained a space and parsed as two. The finding says so rather than leaving
//! the reader to work it out.

use crate::document::{Category, CifValue, CifValueRef, DataBlock, Document};
use crate::lexer::{LexError, Lexer, Quoting, Spanned, Token};
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::io::InputBuffer;
use pdbiox_core::span::{ByteSpan, Position};
use std::borrow::Cow;
use std::collections::HashSet;
use std::marker::PhantomData;
use std::sync::Arc;

/// A document and everything that was wrong with the file it came from.
pub type ParseResult = Result<(Document, Vec<Diagnostic>), Vec<Diagnostic>>;

/// One parsed CIF scalar borrowing text directly from the input buffer.
///
/// This hidden integration type lets domain projections consume parser events
/// without building a lossless [`Document`] or allocating text per cell.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CifScalar<'input> {
    /// Written `.`.
    Inapplicable,
    /// Written `?`.
    Unknown,
    /// Text borrowing the input buffer.
    Text(&'input str),
    /// An exact whole number.
    Integer(i64),
    /// A floating-point number.
    Float(f64),
}

/// Compact identity of one declared column within the active loop.
///
/// Hot sinks resolve this ordinal once from the loop header instead of
/// comparing category and item strings for every value.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub(crate) struct ColumnId(u32);

impl ColumnId {
    pub(crate) fn from_position(position: usize) -> Option<Self> {
        u32::try_from(position).ok().map(Self)
    }

    /// Zero-based position within the active declaration.
    #[must_use]
    pub(crate) const fn position(self) -> usize {
        self.0 as usize
    }
}

impl<'input> From<CifValueRef<'input>> for CifScalar<'input> {
    fn from(value: CifValueRef<'input>) -> Self {
        match value {
            CifValueRef::Inapplicable => Self::Inapplicable,
            CifValueRef::Unknown => Self::Unknown,
            CifValueRef::Text(value) => Self::Text(value),
            CifValueRef::Integer(value) => Self::Integer(value),
            CifValueRef::Float(value) => Self::Float(value),
        }
    }
}

/// Event consumer for allocation-bounded projections over CIF input.
///
/// Category and item names are allocated once per declaration by the parser;
/// each cell text borrow is valid only for its [`Self::value`] callback. This
/// deliberately prevents a projection from retaining an input or decoder
/// buffer and allows the same sink to consume text and binary CIF streams.
///
/// ```compile_fail
/// use pdbiox_cif::{CifEventSink, CifScalar};
/// use pdbiox_core::span::ByteSpan;
///
/// struct RetainsText<'a>(Option<&'a str>);
///
/// impl CifEventSink for RetainsText<'_> {
///     type Output = ();
///
///     fn block(&mut self, _name: &str) {}
///
///     fn value(&mut self, _category: &str, _item: &str, value: CifScalar<'_>, _span: ByteSpan) {
///         if let CifScalar::Text(text) = value {
///             self.0 = Some(text);
///         }
///     }
///
///     fn finish(self) {}
/// }
/// ```
#[doc(hidden)]
pub trait CifEventSink {
    /// Completed projection.
    type Output;

    /// Starts a data block.
    fn block(&mut self, name: &str);

    /// Whether values from one category need scalar decoding and delivery.
    fn accepts_category(&self, _category: &str) -> bool {
        true
    }

    /// Starts a loop after all of its tags have been collected.
    fn begin_loop(&mut self, _tags: &[(Box<str>, Box<str>)]) {}

    /// Consumes one typed scalar.
    fn value(&mut self, category: &str, item: &str, value: CifScalar<'_>, span: ByteSpan);

    /// Completes one scalar or loop row.
    fn end_row(&mut self) {}

    /// Completes an active loop.
    fn end_loop(&mut self) {}

    /// Completes the projection.
    fn finish(self) -> Self::Output;
}

/// Reads a document.
///
/// # Errors
///
/// Returns the findings that stopped the read: text that is not valid, a value
/// or quote that never closed, or a file with no block header.
pub fn parse(input: &InputBuffer) -> ParseResult {
    parse_into(input, DocumentSink::all())
}

/// Parses CIF into a caller-provided event projection.
///
/// The same lexer, row state machine, and diagnostics as [`parse`] are used;
/// only the destination representation differs.
///
/// # Errors
///
/// Returns the syntax findings that stopped parsing.
#[doc(hidden)]
pub fn parse_events<S>(
    input: &InputBuffer,
    sink: S,
) -> Result<(S::Output, Vec<Diagnostic>), Vec<Diagnostic>>
where
    S: CifEventSink,
{
    parse_into(input, EventAdapter(sink))
}

struct EventAdapter<S>(S);

impl<'input, S> ValueSink<'input> for EventAdapter<S>
where
    S: CifEventSink,
{
    type Output = S::Output;

    fn block(&mut self, name: &str) {
        self.0.block(name);
    }

    fn begin_loop(&mut self, tags: &[(Box<str>, Box<str>)]) {
        self.0.begin_loop(tags);
    }

    fn value(
        &mut self,
        _column: ColumnId,
        category: &str,
        item: &str,
        text: &'input str,
        quoting: Quoting,
        span: ByteSpan,
    ) {
        if !self.0.accepts_category(category) {
            return;
        }
        self.0.value(
            category,
            item,
            CifScalar::from(CifValueRef::parse(text, quoting)),
            span,
        );
    }

    fn end_row(&mut self) {
        self.0.end_row();
    }

    fn end_loop(&mut self) {
        self.0.end_loop();
    }

    fn finish(self) -> Self::Output {
        self.0.finish()
    }
}

pub(crate) trait ValueSink<'input> {
    type Output;

    fn block(&mut self, name: &str);

    fn begin_loop(&mut self, _tags: &[(Box<str>, Box<str>)]) {}

    fn value(
        &mut self,
        column: ColumnId,
        category: &str,
        item: &str,
        text: &'input str,
        quoting: Quoting,
        span: ByteSpan,
    );

    fn end_row(&mut self) {}

    fn end_loop(&mut self) {}

    fn finish(self) -> Self::Output;
}

pub(crate) fn parse_into<'input, S>(
    input: &'input InputBuffer,
    sink: S,
) -> Result<(S::Output, Vec<Diagnostic>), Vec<Diagnostic>>
where
    S: ValueSink<'input>,
{
    let mut lexer = match Lexer::new(input.as_bytes()) {
        Ok(lexer) => lexer,
        Err(error) => return Err(vec![lex_finding(error)]),
    };
    let mut state = ParseState::new(sink);

    loop {
        match lexer.next_token() {
            Ok(Some(spanned)) => state.token(spanned),
            Ok(None) => break,
            Err(error) => {
                state.findings.push(lex_finding(error));
                return Err(state.findings.finish());
            }
        }
    }
    state.finish()
}

pub(crate) struct DocumentSink {
    document: Document,
    texts: TextInterner,
    keep: fn(&str) -> bool,
}

impl DocumentSink {
    fn all() -> Self {
        Self::with_filter(|_| true)
    }

    pub(crate) fn with_filter(keep: fn(&str) -> bool) -> Self {
        Self {
            document: Document::new(),
            texts: TextInterner::default(),
            keep,
        }
    }

    pub(crate) fn block(&mut self, name: &str) {
        self.document.push(DataBlock::new(name));
    }

    pub(crate) fn value(
        &mut self,
        category: &str,
        item: &str,
        text: &str,
        quoting: Quoting,
        span: ByteSpan,
    ) {
        if !(self.keep)(category) {
            return;
        }
        let Some(block) = self.document.last_block_mut() else {
            return;
        };
        block.category_mut(category, span).column_mut(item).push(
            CifValue::parse_with_text(text, quoting, |text| self.texts.intern(text)),
            quoting,
        );
    }

    pub(crate) fn finish(self) -> Document {
        self.document
    }

    pub(crate) const fn document(&self) -> &Document {
        &self.document
    }
}

impl<'input> ValueSink<'input> for DocumentSink {
    type Output = Document;

    fn block(&mut self, name: &str) {
        Self::block(self, name);
    }

    fn value(
        &mut self,
        _column: ColumnId,
        category: &str,
        item: &str,
        text: &'input str,
        quoting: Quoting,
        span: ByteSpan,
    ) {
        Self::value(self, category, item, text, quoting, span);
    }

    fn finish(self) -> Self::Output {
        Self::finish(self)
    }
}

/// Splits an item name into its category and its item.
///
/// A name without a dot belongs to no category, which older files use and which
/// is kept under its own name rather than discarded.
#[must_use]
pub fn split_tag(tag: &str) -> (&str, &str) {
    let body = match tag.strip_prefix('_') {
        Some(body) => body,
        None => tag,
    };
    match body.split_once('.') {
        Some((category, item)) => (category, item),
        None => (body, body),
    }
}

include!("document/state.rs");

#[derive(Default)]
struct TextInterner {
    values: HashSet<Arc<str>>,
}

impl TextInterner {
    fn intern(&mut self, text: &str) -> Arc<str> {
        if let Some(value) = self.values.get(text) {
            return Arc::clone(value);
        }
        let value: Arc<str> = text.into();
        self.values.insert(Arc::clone(&value));
        value
    }
}

/// Turns a lexer refusal into a finding.
fn lex_finding(error: LexError) -> Diagnostic {
    match error {
        LexError::UnterminatedQuote(at) => Diagnostic::new(Code::E1101).at(ByteSpan::empty(at)),
        LexError::UnterminatedText(at) => Diagnostic::new(Code::E1102).at(ByteSpan::empty(at)),
        LexError::NotText => Diagnostic::new(Code::E1201)
            .with_message("input is not valid text")
            .at(ByteSpan::empty(Position::START)),
        LexError::PositionOverflow(at) => Diagnostic::new(Code::E1201)
            .with_message("input exceeds the supported source-position range")
            .at(ByteSpan::empty(at)),
    }
}

/// The rows of one category, as a cursor a caller can walk.
///
/// Reading a category row by row rather than item by item is what the lowering
/// stage wants, and doing it through one type keeps the absent-item handling in
/// one place.
#[derive(Clone, Copy, Debug)]
pub struct Rows<'a> {
    category: &'a Category,
    row: usize,
}

impl<'a> Rows<'a> {
    /// Walks the rows of a category.
    #[must_use]
    pub const fn new(category: &'a Category) -> Self {
        Self { category, row: 0 }
    }

    /// The category being walked.
    #[must_use]
    pub const fn category(&self) -> &'a Category {
        self.category
    }

    /// The current row number.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }

    /// Moves to the next row, returning false at the end.
    pub fn advance(&mut self) -> bool {
        self.row += 1;
        self.row < self.category.row_count()
    }

    /// One value of the current row, as text.
    #[must_use]
    pub fn text(&self, item: &str) -> Option<&'a str> {
        self.category.value(item, self.row)?.as_str()
    }

    /// One value of the current row, as an identifier.
    ///
    /// Use this wherever the value names something rather than measures it: the
    /// format writes entity numbers and many chain labels as bare numerals, and
    /// reading those only as text loses them.
    #[must_use]
    pub fn identifier(&self, item: &str) -> Option<Cow<'a, str>> {
        self.category.value(item, self.row)?.as_identifier()
    }

    /// One value of the current row, as a whole number.
    #[must_use]
    pub fn integer(&self, item: &str) -> Option<i64> {
        self.category.value(item, self.row)?.as_integer()
    }

    /// One value of the current row, as a number.
    #[must_use]
    pub fn float(&self, item: &str) -> Option<f64> {
        self.category.value(item, self.row)?.as_float()
    }

    /// Whether one value of the current row was recorded at all.
    #[must_use]
    pub fn is_recorded(&self, item: &str) -> bool {
        self.category
            .value(item, self.row)
            .is_some_and(CifValue::is_recorded)
    }
}

#[cfg(test)]
#[path = "document_tests.rs"]
mod tests;
