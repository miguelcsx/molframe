//! Assembling tokens into a document.
//!
//! Two shapes carry data: a tag followed by one value, and a loop header
//! followed by tags and then rows of values. Everything else is structure.
//!
//! A row count that is not a whole multiple of the item count is the single most
//! common way a hand-edited file breaks, almost always because an unquoted value
//! contained a space and parsed as two. The finding says so rather than leaving
//! the reader to work it out.

use crate::document::{Category, CifValue, DataBlock, Document};
use crate::lexer::{LexError, Lexer, Quoting, Spanned, Token};
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::io::InputBuffer;
use pdbiox_core::span::{ByteSpan, Position};
use std::borrow::Cow;

/// A document and everything that was wrong with the file it came from.
pub type ParseResult = Result<(Document, Vec<Diagnostic>), Vec<Diagnostic>>;

/// Reads a document.
///
/// # Errors
///
/// Returns the findings that stopped the read: text that is not valid, a value
/// or quote that never closed, or a file with no block header.
pub fn parse(input: &InputBuffer) -> ParseResult {
    let mut lexer = match Lexer::new(input.as_bytes()) {
        Ok(lexer) => lexer,
        Err(error) => return Err(vec![lex_finding(error)]),
    };
    let mut state = ParseState::default();

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

/// What the parser is in the middle of.
#[derive(Default)]
struct ParseState {
    document: Document,
    findings: Diagnostics,
    /// Item names of the loop being read, in column order.
    loop_tags: Vec<(Box<str>, Box<str>)>,
    /// Where the loop header began.
    loop_span: Option<ByteSpan>,
    /// How many values of the current loop row have been read.
    loop_cursor: usize,
    /// True between `loop_` and its first value.
    collecting_tags: bool,
    /// The item a bare `tag value` pair is waiting on.
    pending_tag: Option<(Box<str>, Box<str>, ByteSpan)>,
    saw_block: bool,
}

impl ParseState {
    fn token(&mut self, spanned: Spanned<'_>) {
        match spanned.token {
            Token::Block(name) => self.block(name),
            Token::FrameStart(_) | Token::FrameEnd => self.end_loop(),
            Token::Loop => self.begin_loop(spanned.span),
            Token::Tag(tag) => self.tag(tag, spanned.span),
            Token::Value(text, quoting) => self.value(text, quoting, spanned.span),
        }
    }

    fn block(&mut self, name: &str) {
        self.end_loop();
        self.saw_block = true;
        self.document.push(DataBlock::new(name));
    }

    fn begin_loop(&mut self, span: ByteSpan) {
        self.end_loop();
        self.collecting_tags = true;
        self.loop_span = Some(span);
    }

    fn tag(&mut self, tag: &str, span: ByteSpan) {
        let (category, item) = split_tag(tag);
        if self.collecting_tags {
            self.loop_tags.push((category.into(), item.into()));
            return;
        }
        self.end_loop();
        self.pending_tag = Some((category.into(), item.into(), span));
    }

    fn value(&mut self, text: &str, quoting: Quoting, span: ByteSpan) {
        if !self.saw_block {
            self.findings.push(Diagnostic::new(Code::E1106).at(span));
            self.saw_block = true;
            self.document.push(DataBlock::new(""));
        }
        self.collecting_tags = false;

        if let Some((category, item, at)) = self.pending_tag.take() {
            store(&mut self.document, &category, &item, text, quoting, at);
            return;
        }
        if self.loop_tags.is_empty() {
            self.findings.push(Diagnostic::new(Code::E1105).at(span));
            return;
        }
        let column = self.loop_cursor % self.loop_tags.len();
        let Some((category, item)) = self.loop_tags.get(column) else {
            return;
        };
        self.loop_cursor += 1;
        store(&mut self.document, category, item, text, quoting, span);
    }

    /// Closes the loop being read, checking that its rows are whole.
    fn end_loop(&mut self) {
        if !self.loop_tags.is_empty() && !self.loop_cursor.is_multiple_of(self.loop_tags.len()) {
            let columns = self.loop_tags.len();
            let short = self.loop_cursor % columns;
            let mut finding = Diagnostic::new(Code::E1103)
                .with_context("items", columns.to_string())
                .with_context("values in the final row", short.to_string());
            if let Some(span) = self.loop_span {
                finding = finding.at(span);
            }
            if let Some((category, _)) = self.loop_tags.first() {
                finding = finding.in_category(category.to_string());
            }
            self.findings.push(finding);
        }
        self.loop_tags.clear();
        self.loop_cursor = 0;
        self.loop_span = None;
        self.collecting_tags = false;
    }

    fn finish(mut self) -> ParseResult {
        self.end_loop();
        if self.pending_tag.is_some() {
            self.findings.push(
                Diagnostic::new(Code::E1105)
                    .with_message("an item name was not followed by a value"),
            );
        }
        if !self.saw_block {
            self.findings.push(Diagnostic::new(Code::E1106));
            return Err(self.findings.finish());
        }
        Ok((self.document, self.findings.finish()))
    }
}

fn store(
    document: &mut Document,
    category: &str,
    item: &str,
    text: &str,
    quoting: Quoting,
    span: ByteSpan,
) {
    let Some(block) = document.last_block_mut() else {
        return;
    };
    block
        .category_mut(category, span)
        .column_mut(item)
        .push(CifValue::parse(text, quoting), quoting);
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
