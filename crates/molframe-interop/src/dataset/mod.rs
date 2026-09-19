//! Lazy manifest-backed structure datasets.

mod filtering;
mod manifest;
mod splitting;

pub use filtering::DatasetFilter;
pub use manifest::{Dataset, DatasetError, LoadError, ManifestEntry};
pub use splitting::{DatasetSplit, DatasetWarning, SplitOptions, SplitRatios, SplitStrategy};
