//! Amino-acid substitution matrices and their lookup.
//!
//! A substitution matrix assigns an integer score to every ordered pair of the
//! twenty standard amino acids together with the ambiguity codes `B` (Asx),
//! `Z` (Glx), `X` (any) and the `*` stop symbol, laid out in the conventional
//! order `ARNDCQEGHILKMFPSTWYVBZX*`. A residue byte is mapped to its row and
//! column through a fixed 256-entry table, so a lookup is `O(1)` and allocates
//! nothing. Input is folded to upper case first, and a byte outside the
//! alphabet is scored as `X`, the "any residue" column, rather than rejected.

use super::profiles::{MatrixError, MatrixIdentity};
use crate::scoring::Score;

/// The number of symbols in the standard amino-acid substitution alphabet.
const ORDER_LEN: usize = 24;

/// The residue bytes in the order the matrix rows and columns follow.
const ORDER: [u8; ORDER_LEN] = *b"ARNDCQEGHILKMFPSTWYVBZX*";

/// The column used for any byte not otherwise in the alphabet: `X`, "any".
const ANY_INDEX: u8 = 22;

/// A dense substitution matrix over the standard amino-acid alphabet.
///
/// The scores are held as a `24 × 24` grid in the [`ORDER`] layout, with a
/// side table turning a residue byte into its index in that grid.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubstitutionMatrix {
    scores: Vec<i32>,
    side: usize,
    lookup: [u8; 256],
    identity: MatrixIdentity,
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
    /// use molframe_seq::matrix::blosum62;
    /// // The distributed matrices are built through this constructor.
    /// let m = blosum62();
    /// assert_eq!(m.get(b'C', b'C'), 9);
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
            identity: MatrixIdentity::custom(),
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
    /// use molframe_seq::matrix::blosum62;
    /// let m = blosum62();
    /// assert_eq!(m.get(b'W', b'W'), 11);
    /// assert_eq!(m.get(b'w', b'W'), 11);
    /// assert_eq!(m.get(b'A', b'W'), -3);
    /// ```
    #[must_use]
    pub fn get(&self, a: u8, b: u8) -> i32 {
        let ai = usize::from(self.lookup[usize::from(a.to_ascii_uppercase())]);
        let bi = usize::from(self.lookup[usize::from(b.to_ascii_uppercase())]);
        self.scores[ai * self.side + bi]
    }

    /// Stable name, version and source of the selected score table.
    #[must_use]
    pub const fn identity(&self) -> &MatrixIdentity {
        &self.identity
    }

    pub(super) fn from_profile(
        alphabet: &[u8],
        scores: Vec<i32>,
        identity: MatrixIdentity,
    ) -> Result<Self, MatrixError> {
        let side = alphabet.len();
        if side == 0 || side > usize::from(u8::MAX) || scores.len() != side * side {
            return Err(MatrixError::MalformedProfile);
        }
        let unknown = alphabet
            .iter()
            .position(|symbol| matches!(symbol, b'X' | b'N'))
            .ok_or(MatrixError::MissingUnknownSymbol)?;
        let unknown = u8::try_from(unknown).map_err(|_| MatrixError::MalformedProfile)?;
        let mut lookup = [unknown; 256];
        for (index, symbol) in alphabet.iter().copied().enumerate() {
            let index = u8::try_from(index).map_err(|_| MatrixError::MalformedProfile)?;
            lookup[usize::from(symbol)] = index;
            lookup[usize::from(symbol.to_ascii_lowercase())] = index;
        }
        Ok(Self {
            scores,
            side,
            lookup,
            identity,
        })
    }
}

impl Score for SubstitutionMatrix {
    fn score(&self, a: u8, b: u8) -> i32 {
        self.get(a, b)
    }
}

/// The BLOSUM62 substitution matrix.
///
/// BLOSUM62 is the default matrix for protein comparison: scores derived from
/// blocks of aligned proteins clustered at 62 % identity. Larger entries on the
/// diagonal mark residues that are conserved (tryptophan at 11, cysteine at 9);
/// off-diagonal entries are positive where a substitution is common and
/// negative where it is rare.
///
/// # Examples
///
/// ```
/// use molframe_seq::matrix::blosum62;
/// let m = blosum62();
/// assert_eq!(m.get(b'A', b'A'), 4);
/// assert_eq!(m.get(b'W', b'W'), 11);
/// assert_eq!(m.get(b'A', b'W'), -3);
/// ```
#[must_use]
pub fn blosum62() -> SubstitutionMatrix {
    #[rustfmt::skip]
    let scores: [[i32; ORDER_LEN]; ORDER_LEN] = [
        [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4],
        [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4],
        [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4],
        [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4],
        [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
        [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4],
        [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
        [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4],
        [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4],
        [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4],
        [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4],
        [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4],
        [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4],
        [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4],
        [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
        [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4],
        [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4],
        [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4],
        [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4],
        [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4],
        [-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4],
        [-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
        [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4],
        [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1],
    ];
    SubstitutionMatrix::new(scores)
}

#[cfg(test)]
#[path = "table_tests.rs"]
mod tests;
