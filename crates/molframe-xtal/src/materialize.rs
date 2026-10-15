//! Explicit expansion of a lazy biological assembly.

use molframe_core::ChainIndex;
use molframe_geom::Rigid;
use std::collections::BTreeMap;
use std::ops::Range;

/// Describes one materialized copy of a source-chain atom span.
#[derive(Clone, Debug)]
struct CopySpan {
    source_chain: ChainIndex,
    transform: Rigid,
    source_atoms: Range<u32>,
    output_start: u32,
}

/// Groups copy indices by source chain for one identical rigid transform.
type CopyGroup = BTreeMap<u32, Vec<usize>>;

/// Exact bit representation of a rigid transform used for deterministic grouping.
type TransformKey = [u64; 12];

mod copy;
mod errors;
mod expand;
mod labels;
mod stream;

use errors::{capacity, single};
use expand::Materializer;
use labels::LabelLedger;
pub use stream::AssemblyChunks;

#[cfg(test)]
#[path = "materialize_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "materialize_stream_tests.rs"]
mod stream_tests;
