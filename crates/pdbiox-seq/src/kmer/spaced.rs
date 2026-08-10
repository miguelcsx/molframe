//! Explicit spaced-seed masks.

/// A spaced-seed mask has no selected positions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SeedPatternError;

impl std::fmt::Display for SeedPatternError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("spaced-seed pattern must select at least one position")
    }
}

impl std::error::Error for SeedPatternError {}

/// A caller-supplied pattern selecting positions from each source window.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedPattern {
    selected: Vec<usize>,
    span: usize,
}

impl SeedPattern {
    /// Builds a pattern from `true` selected and `false` ignored positions.
    ///
    /// # Errors
    ///
    /// Returns [`SeedPatternError`] when no position is selected.
    pub fn new(mask: &[bool]) -> Result<Self, SeedPatternError> {
        let selected = mask
            .iter()
            .enumerate()
            .filter_map(|(index, selected)| selected.then_some(index))
            .collect::<Vec<_>>();
        if selected.is_empty() {
            Err(SeedPatternError)
        } else {
            Ok(Self {
                selected,
                span: mask.len(),
            })
        }
    }

    /// Width of the source window.
    #[must_use]
    pub const fn span(&self) -> usize {
        self.span
    }

    /// Number of symbols retained in a seed.
    #[must_use]
    pub fn weight(&self) -> usize {
        self.selected.len()
    }

    pub(super) fn seeds<'a>(
        &'a self,
        sequence: &'a [u8],
    ) -> impl Iterator<Item = (usize, Vec<u8>)> + 'a {
        sequence
            .windows(self.span)
            .enumerate()
            .map(|(position, window)| {
                let seed = self.selected.iter().map(|index| window[*index]).collect();
                (position, seed)
            })
    }
}

#[cfg(test)]
#[path = "spaced_tests.rs"]
mod tests;
