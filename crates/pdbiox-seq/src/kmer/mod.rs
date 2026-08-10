//! Exact, spaced, bucketed, and similarity-aware sequence words.

mod counts;
mod similar;
mod spaced;
mod table;

pub use counts::*;
pub use similar::*;
pub use spaced::*;
pub use table::*;
