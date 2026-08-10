//! How much an alignment column is worth.
//!
//! A scoring scheme rewards a match, penalises a mismatch, and charges a gap in
//! two parts: a one-off cost to open it and a smaller cost for each residue it
//! spans. Penalties are held as the values added to the score, so they are
//! zero or negative and a length-`k` gap costs `open + (k - 1) * extend`.

/// A match/mismatch scoring scheme with affine gap costs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scoring {
    /// Score added when two residues are equal.
    pub match_score: i32,
    /// Score added when two residues differ.
    pub mismatch_score: i32,
    /// Score added for the first residue of a gap. Zero or negative.
    pub gap_open: i32,
    /// Score added for each residue of a gap after the first. Zero or negative.
    pub gap_extend: i32,
}

impl Scoring {
    /// A conventional nucleotide-style scheme: +1 match, −1 mismatch, and a
    /// gap that costs −2 to open and −1 to extend.
    #[must_use]
    pub const fn simple() -> Self {
        Self {
            match_score: 1,
            mismatch_score: -1,
            gap_open: -2,
            gap_extend: -1,
        }
    }

    /// The substitution score for a pair of residues.
    #[must_use]
    pub const fn substitution(&self, left: u8, right: u8) -> i32 {
        if left == right {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// Something that can score one residue byte against another.
///
/// This is the single point of contact between an alignment and how a column is
/// valued, so the same dynamic program can run against a plain match/mismatch
/// [`Scoring`] or a full substitution matrix without either knowing about the
/// other. Only the substitution part lives here; gap costs are passed to the
/// aligners separately because they belong to the alignment, not the alphabet.
pub trait Score {
    /// The score contributed by aligning residue `a` against residue `b`.
    fn score(&self, a: u8, b: u8) -> i32;
}

impl Score for Scoring {
    fn score(&self, a: u8, b: u8) -> i32 {
        self.substitution(a, b)
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
