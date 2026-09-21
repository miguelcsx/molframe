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

use memchr::{memchr, memchr_iter, memrchr};
use molframe_core::span::{ByteSpan, Position};

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
/// use molframe_cif::lexer::{Lexer, Token};
///
/// let mut lexer = Lexer::new(b"data_test\n_entry.id 1ABC\n")?;
/// let first = lexer.next_token()?;
/// assert_eq!(first.map(|spanned| spanned.token), Some(Token::Block("test")));
/// # Ok::<(), molframe_cif::lexer::LexError>(())
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

        // Nearly every bare word in a coordinate loop is a value, and the three
        // keywords begin with three distinct letters. Deciding on the first
        // byte turns three prefix comparisons per value into one.
        match word.as_bytes().first() {
            Some(b'd' | b'D') => {
                if let Some(name) = strip_prefix_ignore_case(word, "data_") {
                    return Ok(Token::Block(name));
                }
            }
            Some(b's' | b'S') => {
                if let Some(name) = strip_prefix_ignore_case(word, "save_") {
                    return Ok(if name.is_empty() {
                        Token::FrameEnd
                    } else {
                        Token::FrameStart(name)
                    });
                }
            }
            Some(b'l' | b'L') if word.eq_ignore_ascii_case("loop_") => {
                return Ok(Token::Loop);
            }
            _ => {}
        }
        Ok(Token::Value(word, Quoting::Bare))
    }

    fn peek(&self) -> Option<u8> {
        self.remaining().as_bytes().first().copied()
    }

    fn remaining(&self) -> &'a str {
        match usize::try_from(self.at.byte_offset) {
            Ok(offset) => slice_from(self.input, offset),
            Err(_) => "",
        }
    }

    /// Advances over `count` bytes, keeping the line and column right.
    ///
    /// A long run has its newlines found with a vectorised search; a short one
    /// walks its bytes. The threshold is not a micro-optimisation but the whole
    /// point: most advances in an `atom_site` loop are a single separating
    /// space, and setting up a vectorised search over one byte costs more than
    /// inspecting it.
    fn advance(&mut self, count: usize) -> Result<(), LexError> {
        /// Shortest run for which a vectorised newline search pays for itself.
        const BULK_THRESHOLD: usize = 32;

        let text = self.remaining();
        let skipped = match text.as_bytes().get(..count) {
            Some(skipped) => skipped,
            None => text.as_bytes(),
        };
        if skipped.len() < BULK_THRESHOLD {
            return self.advance_bytewise(skipped);
        }
        self.advance_all(skipped)
    }

    /// Advances over a short run one byte at a time.
    fn advance_bytewise(&mut self, skipped: &[u8]) -> Result<(), LexError> {
        for byte in skipped {
            self.at = self
                .at
                .advance(*byte)
                .ok_or(LexError::PositionOverflow(self.at))?;
        }
        Ok(())
    }

    /// Advances over exactly `skipped`, updating the position in bulk.
    fn advance_all(&mut self, skipped: &[u8]) -> Result<(), LexError> {
        let Ok(length) = u64::try_from(skipped.len()) else {
            return Err(LexError::PositionOverflow(self.at));
        };
        let Some(byte_offset) = self.at.byte_offset.checked_add(length) else {
            return Err(LexError::PositionOverflow(self.at));
        };

        let newlines = memchr_iter(b'\n', skipped).count();
        if newlines == 0 {
            let Ok(narrow) = u64::try_from(skipped.len()) else {
                return Err(LexError::PositionOverflow(self.at));
            };
            let Some(column) = self.at.column.checked_add(narrow) else {
                return Err(LexError::PositionOverflow(self.at));
            };
            self.at.byte_offset = byte_offset;
            self.at.column = column;
            return Ok(());
        }

        let Ok(newlines) = u64::try_from(newlines) else {
            return Err(LexError::PositionOverflow(self.at));
        };
        let Some(line) = self.at.line.checked_add(newlines) else {
            return Err(LexError::PositionOverflow(self.at));
        };
        let Some(last) = memrchr(b'\n', skipped) else {
            return Err(LexError::PositionOverflow(self.at));
        };
        let Ok(after) = u64::try_from(skipped.len() - last - 1) else {
            return Err(LexError::PositionOverflow(self.at));
        };
        let Some(column) = after.checked_add(1) else {
            return Err(LexError::PositionOverflow(self.at));
        };

        self.at.byte_offset = byte_offset;
        self.at.line = line;
        self.at.column = column;
        Ok(())
    }

    /// Advances a run the caller has already established holds no newline.
    ///
    /// Bare values are scanned once to find their end, and that scan proves the
    /// absence of a newline. Re-deriving it here would walk the same bytes a
    /// second time, so this is pure arithmetic — which matters because in an
    /// `atom_site` loop nearly every token takes this path.
    fn advance_non_newline(&mut self, count: usize) -> Result<(), LexError> {
        let Ok(narrow) = u64::try_from(count) else {
            return self.advance(count);
        };
        let Some(byte_offset) = self.at.byte_offset.checked_add(narrow) else {
            return self.advance(count);
        };
        let Some(column) = self.at.column.checked_add(narrow) else {
            return self.advance(count);
        };
        self.at.byte_offset = byte_offset;
        self.at.column = column;
        Ok(())
    }

    /// Skips whitespace and comments.
    ///
    /// A comment runs to the end of its line, and a `#` inside a quoted value is
    /// an ordinary character — which is why this only runs between tokens.
    fn skip_trivia(&mut self) -> Result<(), LexError> {
        loop {
            let text = self.remaining();
            let bytes = text.as_bytes();
            let spaces = match bytes.iter().position(|byte| !byte.is_ascii_whitespace()) {
                Some(spaces) => spaces,
                None => bytes.len(),
            };
            // Advance over the run already located, rather than counting it
            // here and letting `advance` walk the same bytes a second time.
            match bytes.get(..spaces) {
                Some(run) => self.advance_all(run)?,
                None => self.advance(spaces)?,
            }

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
        self.advance_non_newline(end)?;
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
    let Some(slice) = text.get(start..) else {
        return "";
    };
    slice
}

/// The text up to `end`, or nothing when `end` does not fall on a boundary.
fn slice_to(text: &str, end: usize) -> &str {
    let Some(slice) = text.get(..end) else {
        return "";
    };
    slice
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
