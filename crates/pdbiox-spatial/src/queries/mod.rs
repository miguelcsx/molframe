//! Public fixed-radius query operations.

mod count;
mod indices;
mod reduce;
pub(crate) mod search;
mod visit;
mod within;

pub use count::count_pairs_within;
pub use reduce::{PairQuery, reduce_pairs_within_unsorted};
pub use search::{
    pairs_within, pairs_within_unsorted, pairs_within_unsorted_with_options,
    pairs_within_with_options,
};
pub use visit::for_each_pairs_within_unsorted;
pub use within::{within, within_with_options};
