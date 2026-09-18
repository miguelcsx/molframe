//! Public pairwise-alignment values and region policy.

use std::ops::Range;

/// One column of an alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Column {
    /// Index into the first sequence, or `None` where it has a gap.
    pub left: Option<usize>,
    /// Index into the second sequence, or `None` where it has a gap.
    pub right: Option<usize>,
}

/// A completed alignment: its score and the columns that achieve it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    /// The optimal score under the scoring scheme.
    pub score: i32,
    /// The aligned columns, in order from the start of the alignment.
    pub columns: Vec<Column>,
}

/// Why an exact pairwise alignment could not be represented or completed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignError {
    /// A dynamic-programming dimension or score exceeded its exact numeric domain.
    NumericOverflow,
    /// Explicit constraints exclude every complete alignment path.
    NoAlignmentPath,
}

impl std::fmt::Display for AlignError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NumericOverflow => formatter
                .write_str("alignment dimensions or score exceed the supported exact range"),
            Self::NoAlignmentPath => {
                formatter.write_str("alignment constraints exclude every complete path")
            }
        }
    }
}

impl std::error::Error for AlignError {}

/// A region is invalid or its exact alignment score is not representable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RegionAlignError {
    /// One requested source range is invalid.
    Region(RegionError),
    /// The exact alignment exceeded the supported numeric domain.
    Alignment(AlignError),
}

impl std::fmt::Display for RegionAlignError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Region(error) => std::fmt::Display::fmt(error, formatter),
            Self::Alignment(error) => std::fmt::Display::fmt(error, formatter),
        }
    }
}

impl std::error::Error for RegionAlignError {}

impl From<RegionError> for RegionAlignError {
    fn from(error: RegionError) -> Self {
        Self::Region(error)
    }
}

impl From<AlignError> for RegionAlignError {
    fn from(error: AlignError) -> Self {
        Self::Alignment(error)
    }
}

/// Endpoint policy for a region-restricted pairwise alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignmentMode {
    /// Both selected regions are aligned end to end.
    Global,
    /// The best-scoring subalignment inside both regions is returned.
    Local,
    /// Terminal gaps inside the selected regions are free.
    SemiGlobal,
}

/// Explicit ranges and search mode for region-restricted alignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionOptions {
    /// Half-open range in the left sequence.
    pub left: Range<usize>,
    /// Half-open range in the right sequence.
    pub right: Range<usize>,
    /// Endpoint policy within those ranges.
    pub mode: AlignmentMode,
    /// Optional diagonal width relative to the selected region starts.
    pub band: Option<usize>,
}

/// A selected sequence range lies outside its source sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RegionError {
    /// Whether the invalid range belongs to the left sequence.
    pub left: bool,
    /// Source sequence length.
    pub length: usize,
}

impl std::fmt::Display for RegionError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let side = if self.left { "left" } else { "right" };
        write!(
            formatter,
            "{side} alignment region exceeds sequence length {}",
            self.length
        )
    }
}

impl std::error::Error for RegionError {}
