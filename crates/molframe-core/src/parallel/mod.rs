//! Dividing work so that the answer does not depend on how many threads ran it.
//!
//! This crate owns the vocabulary and, behind the off-by-default `rayon`
//! feature, the executor. It does not parallelise its own operations: nothing
//! here runs a kernel, and nothing above the type layer decides a worker count.
//! That decision belongs to the caller (NFR-113).

mod consume;
mod execute;
mod ordered;
mod plan;
mod policy;

pub use consume::{BlockExecutionError, try_for_each_block_in};
pub(crate) use execute::SharedPool;
pub use execute::{WorkerPanicked, map_blocks_in};
pub use plan::{BlockPlan, DEFAULT_BLOCK_ITEMS};
pub use policy::ReductionPolicy;
