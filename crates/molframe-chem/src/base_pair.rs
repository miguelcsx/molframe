//! Canonical nucleotide identity, edge contacts and reference geometry.

mod edges;
pub use edges::{WatsonCrickGeometry, has_watson_crick_geometry, watson_crick_contact};

/// Canonical identity used by edge contacts and dictionary-backed analysis.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CanonicalBase {
    /// Adenine.
    Adenine,
    /// Cytosine.
    Cytosine,
    /// Guanine.
    Guanine,
    /// Thymine.
    Thymine,
    /// Uracil.
    Uracil,
}

impl CanonicalBase {
    /// Resolves a canonical one-letter code.
    #[must_use]
    pub const fn from_one_letter(code: u8) -> Option<Self> {
        match code.to_ascii_uppercase() {
            b'A' => Some(Self::Adenine),
            b'C' => Some(Self::Cytosine),
            b'G' => Some(Self::Guanine),
            b'T' => Some(Self::Thymine),
            b'U' => Some(Self::Uracil),
            _ => None,
        }
    }

    /// Whether two bases form a canonical complementary pair.
    #[must_use]
    pub const fn complementary(self, other: Self) -> bool {
        matches!(
            (self, other),
            (Self::Adenine, Self::Thymine | Self::Uracil)
                | (Self::Thymine | Self::Uracil, Self::Adenine)
                | (Self::Guanine, Self::Cytosine)
                | (Self::Cytosine, Self::Guanine)
        )
    }
}
