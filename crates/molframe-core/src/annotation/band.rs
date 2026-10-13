//! Named pLDDT confidence bands shared by selection, colouring and bindings.
//!
//! The four bands are the ones structure viewers report for predicted models.
//! Boundaries are half-open from below so every finite score belongs to exactly
//! one band: `very_high` is `> 90`, `confident` is `(70, 90]`, `low` is
//! `(50, 70]` and `very_low` is `<= 50`.

/// One of the four pLDDT confidence bands.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum PlddtBand {
    /// Score of at most 50.
    VeryLow,
    /// Score above 50 and at most 70.
    Low,
    /// Score above 70 and at most 90.
    Confident,
    /// Score above 90.
    VeryHigh,
}

impl PlddtBand {
    /// All bands from least to most confident.
    pub const ALL: [Self; 4] = [Self::VeryLow, Self::Low, Self::Confident, Self::VeryHigh];

    /// Band of `score`, or `None` for a non-finite score.
    #[must_use]
    pub fn classify(score: f64) -> Option<Self> {
        if !score.is_finite() {
            return None;
        }
        Some(if score > 90.0 {
            Self::VeryHigh
        } else if score > 70.0 {
            Self::Confident
        } else if score > 50.0 {
            Self::Low
        } else {
            Self::VeryLow
        })
    }

    /// Stable lowercase name, also the selection macro suffix.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::VeryLow => "very_low",
            Self::Low => "low",
            Self::Confident => "confident",
            Self::VeryHigh => "very_high",
        }
    }

    /// Parses a name produced by [`Self::name`], ignoring ASCII case.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|band| band.name().eq_ignore_ascii_case(name))
    }

    /// Conventional display colour as sRGB bytes.
    #[must_use]
    pub const fn rgb(self) -> [u8; 3] {
        match self {
            Self::VeryLow => [255, 125, 69],
            Self::Low => [255, 219, 19],
            Self::Confident => [101, 203, 243],
            Self::VeryHigh => [0, 83, 214],
        }
    }

    /// Lowest score that can belong to this band, exclusive except for `VeryLow`.
    #[must_use]
    pub const fn lower_bound(self) -> f64 {
        match self {
            Self::VeryLow => f64::NEG_INFINITY,
            Self::Low => 50.0,
            Self::Confident => 70.0,
            Self::VeryHigh => 90.0,
        }
    }
}

#[cfg(test)]
#[path = "band_tests.rs"]
mod tests;
