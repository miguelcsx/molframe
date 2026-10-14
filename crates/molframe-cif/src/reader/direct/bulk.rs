//! Speculative parallel reading of an `atom_site` loop body.
//!
//! A CIF token never spans a line break except inside a `;` text block, so the
//! body can be cut at line starts and each piece lexed on its own. The cuts are
//! speculation: a piece is trusted only when every piece before it was whole
//! rows, and the first piece that is not (a text block, a lex error, a count
//! that is not a multiple of the column count) ends the bulk read there. The
//! ordinary one-value-at-a-time parser resumes from that point, so findings,
//! spans and row order are exactly those of a serial read.

use super::row::AtomRow;
use crate::document::CifValueRef;
use crate::lexer::{Lexer, Quoting, Token};
use crate::lower::Field;
use molframe_core::execution::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};

/// Bytes handed to the workers at a time, bounding the rows held at once.
const WINDOW_BYTES: usize = 4 * 1024 * 1024;
/// Bytes one worker task lexes.
const CHUNK_BYTES: usize = 64 * 1024;

struct Chunk<'input> {
    rows: Vec<AtomRow<'input>>,
    /// Offset of the first token that is not a value, when the piece has one.
    terminator: Option<usize>,
    whole: bool,
}

impl Chunk<'_> {
    const fn bad() -> Self {
        Self {
            rows: Vec::new(),
            terminator: None,
            whole: false,
        }
    }
}

fn scan<'input>(text: &'input str, columns: &[Option<Field>]) -> Chunk<'input> {
    let width = columns.len();
    if width == 0 {
        return Chunk::bad();
    }
    let mut lexer = Lexer::over(text);
    let mut rows = Vec::with_capacity(text.len() / 48);
    let mut row = AtomRow::new();
    let mut cursor = 0;
    let mut terminator = None;
    loop {
        let spanned = match lexer.next_token() {
            Ok(Some(spanned)) => spanned,
            Ok(None) => break,
            Err(_) => return Chunk::bad(),
        };
        match spanned.token {
            Token::Value(_, Quoting::Text) => return Chunk::bad(),
            Token::Value(value, quoting) => {
                if let Some(Some(field)) = columns.get(cursor) {
                    row.set(*field, CifValueRef::parse(value, quoting));
                }
                cursor += 1;
                if cursor == width {
                    cursor = 0;
                    rows.push(std::mem::replace(&mut row, AtomRow::new()));
                }
            }
            _ => {
                let Ok(start) = usize::try_from(spanned.span.start.byte_offset) else {
                    return Chunk::bad();
                };
                terminator = Some(start);
                break;
            }
        }
    }
    Chunk {
        rows,
        terminator,
        whole: cursor == 0,
    }
}

/// Start offsets of the pieces of `window`, ending with its length.
fn cuts(window: &str) -> Vec<usize> {
    let bytes = window.as_bytes();
    let mut cuts = vec![0];
    let mut at: usize = 0;
    while let Some(next) = at
        .checked_add(CHUNK_BYTES)
        .filter(|next| *next < bytes.len())
    {
        let Some(newline) = bytes
            .get(next..)
            .and_then(|tail| memchr::memchr(b'\n', tail))
        else {
            break;
        };
        at = next + newline + 1;
        if at >= bytes.len() {
            break;
        }
        cuts.push(at);
    }
    cuts.push(bytes.len());
    cuts
}

/// Hands whole leading rows of `text` to `apply` in order and returns how many
/// bytes they cover, which is always a row boundary.
pub(super) fn run<'input>(
    text: &'input str,
    columns: &[Option<Field>],
    mut apply: impl FnMut(AtomRow<'input>),
) -> usize {
    let context = ExecutionContext::default();
    let mut done = 0;
    while done < text.len() {
        let mut end = done.saturating_add(WINDOW_BYTES).min(text.len());
        if end < text.len() {
            let Some(newline) = text
                .as_bytes()
                .get(end..)
                .and_then(|tail| memchr::memchr(b'\n', tail))
            else {
                return done;
            };
            end += newline + 1;
        }
        let Some(window) = text.get(done..end) else {
            return done;
        };
        let cuts = cuts(window);
        let plan = BlockPlan::new(cuts.len() - 1, 1);
        let Ok(chunks) = map_blocks_in(plan, &context, |index, _| {
            let range = cuts.get(index).copied().zip(cuts.get(index + 1).copied());
            let piece = range.and_then(|(from, to)| window.get(from..to));
            match piece {
                Some(piece) => scan(piece, columns),
                None => Chunk::bad(),
            }
        }) else {
            return done;
        };
        for (chunk, start) in chunks.into_iter().zip(cuts.iter().copied()) {
            if !chunk.whole {
                return done + start;
            }
            for row in chunk.rows {
                apply(row);
            }
            if let Some(offset) = chunk.terminator {
                return done + start + offset;
            }
        }
        done = end;
    }
    done
}
