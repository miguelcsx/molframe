//! A residue substitution matrix as a reference table: alphabet row, then one row per symbol.
//!
//! The scores are published tables (BLOSUM, PAM, NUC.4.4, ...), so they are files
//! selected by name like every other reference table (profile line `score-matrix`,
//! files `<name>.matrix`), with their identity recorded in results; nothing is
//! compiled in. A file lists the alphabet once, then a labelled row of integers
//! per symbol:
//!
//! ```text
//! #  A  R  N
//! A  4 -1 -2
//! R -1  5  0
//! N -2  0  6
//! ```
//!
//! The alphabet must contain an `X` or `N` ambiguity symbol: a byte outside the
//! alphabet is scored as that symbol rather than refused.
//!
//! Cost: decoding is `O(side^2)` once per table; a lookup is two array reads.

use molframe_chem::{
    Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};

const WHAT: &str = "score matrix";

/// A validated square score table over a byte alphabet.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScoreMatrix {
    pub(super) alphabet: Vec<u8>,
    pub(super) scores: Vec<i32>,
    rows: usize,
}

impl TableKind for ScoreMatrix {
    const KIND: &'static str = "score-matrix";
    const EXTENSION: &'static str = "matrix";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(WHAT, text_of(WHAT, bytes)?, &mut table)?;
        if table.alphabet.is_empty() {
            return Err(TableTextError::at(WHAT, 0, "the table lists no alphabet"));
        }
        if table.rows != table.alphabet.len() {
            return Err(TableTextError::at(
                WHAT,
                0,
                format!(
                    "{} symbols but {} rows; the matrix must be square",
                    table.alphabet.len(),
                    table.rows
                ),
            ));
        }
        if !table
            .alphabet
            .iter()
            .any(|symbol| matches!(symbol, b'X' | b'N'))
        {
            return Err(TableTextError::at(
                WHAT,
                0,
                "the alphabet has no `X` or `N` ambiguity symbol",
            ));
        }
        Ok(Provided { identity, table })
    }
}

impl TableLines for ScoreMatrix {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let row = Row::new(WHAT, line);
        let fail = |message: String| row.fail(message);
        if self.alphabet.is_empty() {
            for word in text.split_whitespace() {
                let &[symbol] = word.as_bytes() else {
                    return Err(fail(format!("`{word}` is not a single-character symbol")));
                };
                if self.alphabet.contains(&symbol) {
                    return Err(fail(format!("symbol `{word}` is listed twice")));
                }
                self.alphabet.push(symbol);
            }
            if self.alphabet.len() > usize::from(u8::MAX) {
                return Err(fail("the alphabet is longer than 255 symbols".into()));
            }
            return Ok(());
        }
        let mut words = text.split_whitespace();
        let label = words.next().map(str::as_bytes);
        if label != self.alphabet.get(self.rows).map(std::slice::from_ref) {
            return Err(fail(format!(
                "row {} must be labelled with symbol `{}`",
                self.rows + 1,
                self.alphabet
                    .get(self.rows)
                    .map_or('?', |symbol| char::from(*symbol))
            )));
        }
        let before = self.scores.len();
        for word in words {
            let score = word
                .parse::<i32>()
                .map_err(|_| fail(format!("`{word}` is not an integer score")))?;
            self.scores.push(score);
        }
        if self.scores.len() - before != self.alphabet.len() {
            return Err(fail(format!(
                "expected {} scores, found {}",
                self.alphabet.len(),
                self.scores.len() - before
            )));
        }
        self.rows += 1;
        Ok(())
    }
}
