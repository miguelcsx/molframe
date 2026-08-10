//! Turning bytes into tokens.
//!
//! Hand-written over `&[u8]` rather than assembled from combinators, for one
//! reason that matters more than speed: the cursor already knows its offset,
//! line and column, so every token carries an exact location for free. A parser
//! built on top of it can point at the character that went wrong instead of
//! naming the file.
//!
//! Nothing here allocates. A token borrows its text from the input, and the
//! quoting style travels alongside so that a later stage can tell a value that
//! was written `.` from one that was written `'.'`.

use memchr::memchr;
use pdbiox_core::span::{ByteSpan, Position};

/// How a value was written, which decides what it means.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Quoting {
    /// Written without quotes.
    Bare,
    /// Written between single quotes.
    Single,
    /// Written between double quotes.
    Double,
    /// Written as a semicolon-delimited block.
    Text,
}

/// What the lexer found.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Token<'a> {
    /// The start of a data block, carrying its name.
    Block(&'a str),
    /// The start of a save frame, carrying its name.
    FrameStart(&'a str),
    /// The end of a save frame.
    FrameEnd,
    /// A loop header.
    Loop,
    /// An item name, including its leading underscore.
    Tag(&'a str),
    /// A value, with how it was written.
    Value(&'a str, Quoting),
}

/// A token and where it came from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Spanned<'a> {
    /// What was found.
    pub token: Token<'a>,
    /// Where it was found.
    pub span: ByteSpan,
}

/// Something the lexer could not get past.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LexError {
    /// A quoted value ran to the end of its line without closing.
    UnterminatedQuote(Position),
    /// A semicolon-delimited block ran to the end of the file without closing.
    UnterminatedText(Position),
    /// The input is not valid text.
    NotText,
    /// A source location exceeded the core's 32-bit position representation.
    PositionOverflow(Position),
}

/// Walks the tokens of a file.
///
/// # Examples
///
/// ```
/// use pdbiox_cif::lexer::{Lexer, Token};
///
/// let mut lexer = Lexer::new(b"data_test\n_entry.id 1ABC\n")?;
/// let first = lexer.next_token()?;
/// assert_eq!(first.map(|spanned| spanned.token), Some(Token::Block("test")));
/// # Ok::<(), pdbiox_cif::lexer::LexError>(())
/// ```
#[derive(Clone, Debug)]
pub struct Lexer<'a> {
    input: &'a str,
    at: Position,
}

impl<'a> Lexer<'a> {
    /// Creates a lexer over `bytes`.
    ///
    /// # Errors
    ///
    /// Returns [`LexError::NotText`] when the bytes are not valid text.
    pub fn new(bytes: &'a [u8]) -> Result<Self, LexError> {
        let input = str::from_utf8(bytes).map_err(|_| LexError::NotText)?;
        Ok(Self {
            input,
            at: Position::START,
        })
    }

    /// Where the lexer has reached.
    #[must_use]
    pub const fn position(&self) -> Position {
        self.at
    }

    /// Reads the next token, or `None` at the end of the input.
    ///
    /// # Errors
    ///
    /// Returns the reason the input could not be got past.
    pub fn next_token(&mut self) -> Result<Option<Spanned<'a>>, LexError> {
        self.skip_trivia()?;
        let start = self.at;
        let Some(byte) = self.peek() else {
            return Ok(None);
        };

        let token = match byte {
            b'_' => Token::Tag(self.take_bare()?),
            b'\'' => Token::Value(self.take_quoted(b'\'')?, Quoting::Single),
            b'"' => Token::Value(self.take_quoted(b'"')?, Quoting::Double),
            b';' if start.column == 1 => Token::Value(self.take_text()?, Quoting::Text),
            _ => self.keyword_or_value()?,
        };
        Ok(Some(Spanned {
            token,
            span: ByteSpan::new(start, self.at.byte_offset),
        }))
    }

    /// A bare word, which may be a keyword or an ordinary value.
    fn keyword_or_value(&mut self) -> Result<Token<'a>, LexError> {
        let word = self.take_bare()?;
        if let Some(name) = strip_prefix_ignore_case(word, "data_") {
            return Ok(Token::Block(name));
        }
        if let Some(name) = strip_prefix_ignore_case(word, "save_") {
            return Ok(if name.is_empty() {
                Token::FrameEnd
            } else {
                Token::FrameStart(name)
            });
        }
        if word.eq_ignore_ascii_case("loop_") {
            return Ok(Token::Loop);
        }
        Ok(Token::Value(word, Quoting::Bare))
    }

    fn peek(&self) -> Option<u8> {
        self.remaining().as_bytes().first().copied()
    }

    fn remaining(&self) -> &'a str {
        slice_from(self.input, self.at.byte_offset as usize)
    }

    /// Advances over `count` bytes, keeping the line and column right.
    fn advance(&mut self, count: usize) -> Result<(), LexError> {
        let text = self.remaining();
        for byte in text.as_bytes().iter().take(count) {
            self.at = self
                .at
                .advance(*byte)
                .ok_or(LexError::PositionOverflow(self.at))?;
        }
        Ok(())
    }

    /// Skips whitespace and comments.
    ///
    /// A comment runs to the end of its line, and a `#` inside a quoted value is
    /// an ordinary character — which is why this only runs between tokens.
    fn skip_trivia(&mut self) -> Result<(), LexError> {
        loop {
            let text = self.remaining();
            let spaces = text.bytes().take_while(u8::is_ascii_whitespace).count();
            self.advance(spaces)?;

            if self.peek() != Some(b'#') {
                return Ok(());
            }
            let rest = self.remaining();
            let line = match memchr(b'\n', rest.as_bytes()) {
                Some(line) => line,
                None => rest.len(),
            };
            self.advance(line)?;
        }
    }

    /// A run of non-whitespace.
    fn take_bare(&mut self) -> Result<&'a str, LexError> {
        let text = self.remaining();
        let end = text
            .bytes()
            .take_while(|byte| !byte.is_ascii_whitespace())
            .count();
        self.advance(end)?;
        Ok(slice_to(text, end))
    }

    /// A value between quotes.
    ///
    /// The closing quote is the one followed by whitespace or the end of the
    /// line, so an apostrophe inside a quoted word does not end it early.
    fn take_quoted(&mut self, quote: u8) -> Result<&'a str, LexError> {
        let opened_at = self.at;
        self.advance(1)?;
        let text = self.remaining();
        let bytes = text.as_bytes();

        let mut cursor = 0usize;
        while cursor < bytes.len() {
            let Some(offset) = memchr(quote, bytes_from(bytes, cursor)) else {
                break;
            };
            let candidate = cursor + offset;
            let closes = match bytes.get(candidate + 1) {
                None => true,
                Some(next) => next.is_ascii_whitespace(),
            };
            if closes {
                self.advance(candidate + 1)?;
                return Ok(slice_to(text, candidate));
            }
            if bytes.get(candidate).is_some_and(|byte| *byte == b'\n') {
                break;
            }
            cursor = candidate + 1;
        }
        Err(LexError::UnterminatedQuote(opened_at))
    }

    /// A semicolon-delimited block.
    ///
    /// It begins with a semicolon in the first column and ends with the next
    /// one, so a semicolon anywhere else is ordinary text.
    fn take_text(&mut self) -> Result<&'a str, LexError> {
        let opened_at = self.at;
        self.advance(1)?;
        // The rest of the opening line belongs to the value.
        let text = self.remaining();
        let bytes = text.as_bytes();

        let mut cursor = 0usize;
        loop {
            let Some(offset) = memchr(b'\n', bytes_from(bytes, cursor)) else {
                return Err(LexError::UnterminatedText(opened_at));
            };
            let line_start = cursor + offset + 1;
            if bytes.get(line_start) == Some(&b';') {
                let value = slice_to(text, cursor + offset);
                self.advance(line_start + 1)?;
                return Ok(value.trim_start_matches(['\r', '\n']));
            }
            cursor = line_start;
        }
    }
}

/// The text from `start` onwards, or nothing when `start` is past the end.
fn slice_from(text: &str, start: usize) -> &str {
    match text.get(start..) {
        Some(slice) => slice,
        None => "",
    }
}

/// The text up to `end`, or nothing when `end` does not fall on a boundary.
fn slice_to(text: &str, end: usize) -> &str {
    match text.get(..end) {
        Some(slice) => slice,
        None => "",
    }
}

/// The bytes from `start` onwards, or nothing when `start` is past the end.
fn bytes_from(bytes: &[u8], start: usize) -> &[u8] {
    match bytes.get(start..) {
        Some(slice) => slice,
        None => &[],
    }
}

fn strip_prefix_ignore_case<'a>(text: &'a str, prefix: &str) -> Option<&'a str> {
    let head = text.get(..prefix.len())?;
    head.eq_ignore_ascii_case(prefix)
        .then(|| text.get(prefix.len()..))?
}

#[cfg(test)]
#[path = "tokenize_tests.rs"]
mod tests;
