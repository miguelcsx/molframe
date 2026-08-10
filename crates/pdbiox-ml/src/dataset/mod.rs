//! Lazy manifest-backed structure datasets.

mod filter;
mod manifest;
mod split;

pub use filter::DatasetFilter;
pub use manifest::{Dataset, DatasetError, LoadError, ManifestEntry};
pub use split::{DatasetSplit, DatasetWarning, SplitOptions, SplitRatios, SplitStrategy};

#[cfg(test)]
#[path = "dataset_tests.rs"]
mod tests;
