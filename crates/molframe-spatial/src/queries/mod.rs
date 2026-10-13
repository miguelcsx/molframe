//! Public fixed-radius query operations.

mod count;
mod crowding;
mod indices;
mod knn;
mod knn_between;
mod reduce;
pub(crate) mod search;
mod visit;
mod within;

pub use count::count_pairs_within;
pub use knn::knn_within;
pub use knn_between::knn_between;
pub use reduce::{PairQuery, reduce_pairs_within_unsorted};
pub use search::{
    pairs_within, pairs_within_unsorted, pairs_within_unsorted_with_options,
    pairs_within_with_options,
};
pub use visit::for_each_pairs_within_unsorted;
pub use within::{within, within_with_options};

#[cfg(test)]
mod knn_tests;

#[cfg(test)]
mod knn_between_tests;
