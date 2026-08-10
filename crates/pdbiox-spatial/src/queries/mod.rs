//! Public fixed-radius query operations.

pub(crate) mod search;

pub use search::{pairs_within, pairs_within_with_options, within, within_with_options};
