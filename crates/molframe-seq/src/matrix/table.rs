//! Amino-acid substitution matrices and their lookup.
//!
//! A substitution matrix assigns an integer score to every ordered pair of the
//! twenty standard amino acids together with the ambiguity codes `B` (Asx),
//! `Z` (Glx), `X` (any) and the `*` stop symbol, laid out in the conventional
//! order `ARNDCQEGHILKMFPSTWYVBZX*`. A residue byte is mapped to its row and
//! column through a fixed 256-entry table, so a lookup is `O(1)` and allocates
//! nothing. Input is folded to upper case first, and a byte outside the
//! alphabet is scored as `X`, the "any residue" column, rather than rejected.

use super::ScoreMatrix;
use crate::scoring::Score;
use molframe_chem::{Provided, TableIdentity};

/// The number of symbols in the standard amino-acid substitution alphabet.
const ORDER_LEN: usize = 24;

/// The residue bytes in the order the matrix rows and columns follow.
const ORDER: [u8; ORDER_LEN] = *b"ARNDCQEGHILKMFPSTWYVBZX*";

/// The column used for any byte not otherwise in the alphabet: `X`, "any".
const ANY_INDEX: u8 = 22;

/// A dense substitution matrix over the standard amino-acid alphabet.
///
/// The scores are held as a `24 × 24` grid in the standard order
/// `ARNDCQEGHILKMFPSTWYVBZX*`, with a side table turning a residue byte into
/// its index in that grid.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubstitutionMatrix {
    scores: Vec<i32>,
    side: usize,
    lookup: [u8; 256],
    identity: TableIdentity,
}

impl SubstitutionMatrix {
    /// Builds a matrix from scores given row by row in the standard
    /// `ARNDCQEGHILKMFPSTWYVBZX*` order.
    ///
    /// The byte-to-index table is derived once here so that later lookups are a
    /// single array read.
    ///
    /// # Examples
    ///
    /// ```
    /// use molframe_seq::SubstitutionMatrix;
    /// let mut scores = [[-1; 24]; 24];
    /// for (index, row) in scores.iter_mut().enumerate() {
    ///     row[index] = 5;
    /// }
    /// let m = SubstitutionMatrix::new(scores);
    /// assert_eq!(m.get(b'C', b'C'), 5);
    /// assert_eq!(m.get(b'C', b'W'), -1);
    /// ```
    #[must_use]
    pub fn new(scores: [[i32; ORDER_LEN]; ORDER_LEN]) -> Self {
        let mut lookup = [ANY_INDEX; 256];
        let mut index = 0;
        let mut code = 0u8;
        while index < ORDER_LEN {
            lookup[usize::from(ORDER[index])] = code;
            index += 1;
            code += 1;
        }
        Self {
            scores: scores.into_iter().flatten().collect(),
            side: ORDER_LEN,
            lookup,
            identity: TableIdentity::builtin("custom", "caller-supplied", "caller", "caller"),
        }
    }

    /// The score for aligning residue `a` against residue `b`.
    ///
    /// Lower-case input is folded to upper case, and any byte outside the
    /// alphabet is treated as `X`.
    ///
    /// # Examples
    ///
    /// ```
    /// use molframe_seq::SubstitutionMatrix;
    /// let m = SubstitutionMatrix::new([[2; 24]; 24]);
    /// assert_eq!(m.get(b'w', b'W'), 2);
    /// assert_eq!(m.get(b'?', b'W'), 2);
    /// ```
    #[must_use]
    pub fn get(&self, a: u8, b: u8) -> i32 {
        let ai = usize::from(self.lookup[usize::from(a.to_ascii_uppercase())]);
        let bi = usize::from(self.lookup[usize::from(b.to_ascii_uppercase())]);
        self.scores[ai * self.side + bi]
    }

    /// Stable name, version and source of the selected score table.
    #[must_use]
    pub const fn identity(&self) -> &TableIdentity {
        &self.identity
    }

    /// Builds the dense lookup from a decoded table and keeps its identity.
    #[must_use]
    pub fn from_table(provided: &Provided<ScoreMatrix>) -> Self {
        let table = &provided.table;
        let side = table.alphabet.len();
        let mut lookup = [0u8; 256];
        let mut unknown = 0u8;
        for (code, symbol) in (0u8..=u8::MAX).zip(table.alphabet.iter().copied()) {
            if matches!(symbol, b'X' | b'N') {
                unknown = code;
            }
        }
        lookup.fill(unknown);
        for (code, symbol) in (0u8..=u8::MAX).zip(table.alphabet.iter().copied()) {
            lookup[usize::from(symbol)] = code;
            lookup[usize::from(symbol.to_ascii_lowercase())] = code;
        }
        Self {
            scores: table.scores.clone(),
            side,
            lookup,
            identity: provided.identity.clone(),
        }
    }
}

impl Score for SubstitutionMatrix {
    fn score(&self, a: u8, b: u8) -> i32 {
        self.get(a, b)
    }
}

#[cfg(test)]
#[path = "table_tests.rs"]
mod tests;
