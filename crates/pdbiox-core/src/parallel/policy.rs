//! Whether a reduction may reassociate.
//!
//! Determinism is not free of consequence but it is nearly free of cost: it
//! holds one partial per block live until the merge, and it fixes the order
//! those partials combine in. What it buys is that a result does not change
//! when a caller changes the thread count, which is what makes a recorded
//! result reproducible by someone who has only the record.
//!
//! A caller who would rather have work stealing says so, and the result says
//! so too. A speedup bought by silently changing the answer is not a speedup.

use std::fmt;

/// Whether a parallel reduction must combine partials in block order.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub enum ReductionPolicy {
    /// Partials combine in ascending block index, whatever order they finish in.
    ///
    /// Output is bit-identical at any worker count.
    #[default]
    Deterministic,
    /// Partials combine as they complete.
    ///
    /// Permits work stealing, and reassociates floating-point sums. A result
    /// produced under this policy MUST record it in its provenance.
    Fast,
}

impl ReductionPolicy {
    /// Returns true when combination order is fixed by block index.
    #[must_use]
    pub const fn is_deterministic(self) -> bool {
        matches!(self, Self::Deterministic)
    }

    /// The name recorded in provenance.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Deterministic => "deterministic",
            Self::Fast => "fast",
        }
    }
}

impl fmt::Display for ReductionPolicy {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

#[cfg(test)]
#[path = "policy_tests.rs"]
mod tests;
