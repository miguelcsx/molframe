//! Public fixed-radius query operations.

pub(crate) mod search;

pub use search::{
    for_each_pairs_within_unsorted, pairs_within, pairs_within_unsorted,
    pairs_within_unsorted_with_options, pairs_within_with_options, within, within_with_options,
};
