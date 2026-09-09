//! The contact network shared by every elastic network model.
//!
//! Both the Gaussian and the anisotropic model spring the same graph: sites
//! within a cutoff are connected, and everything that differs between them is
//! what each edge contributes to its operator. Construction is `O(N + E)` for
//! `N` selected sites and `E` cutoff contacts, and every allocation on the path
//! is checked against the caller's ceiling before it is made, so a network that
//! will not fit is refused instead of exhausting the host.

#[path = "network/graph.rs"]
mod graph;

pub(crate) use graph::{ContactGraph, expand_index};

use pdbiox_spatial::{SpatialBackend, SpatialError};

/// What a caller allows the network build to cost, and how it searches.
///
/// Both model families carry richer public option types; this is the subset the
/// shared construction actually reads, so neither model's public surface has to
/// know about the other's.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct NetworkBudget {
    /// Maximum separation defining a network spring.
    pub(crate) contact_distance: f32,
    /// Maximum bytes for graph construction, solver workspace, and output.
    pub(crate) memory_limit_bytes: usize,
    /// Spatial implementation used to build contacts.
    pub(crate) backend: SpatialBackend,
}

/// Why a contact network could not be built within its budget.
#[derive(Debug)]
pub(crate) enum NetworkError {
    /// Graph construction and workspace exceed the caller's ceiling.
    MemoryLimit {
        /// Required workspace bytes.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
    /// The neighbour search failed or an index left the addressable range.
    Spatial(SpatialError),
}

impl From<SpatialError> for NetworkError {
    fn from(error: SpatialError) -> Self {
        Self::Spatial(error)
    }
}

/// Multiplies sizes, reporting an overflow as the ceiling being exceeded.
///
/// A product that cannot be represented is by definition larger than any
/// ceiling a caller could have set, so it is reported the same way rather than
/// as a separate arithmetic failure the caller would have to handle twice.
pub(crate) fn checked_product(
    factors: &[usize],
    budget: NetworkBudget,
) -> Result<usize, NetworkError> {
    factors
        .iter()
        .try_fold(1_usize, |product, factor| product.checked_mul(*factor))
        .ok_or(NetworkError::MemoryLimit {
            required: usize::MAX,
            limit: budget.memory_limit_bytes,
        })
}

/// Adds sizes, reporting an overflow as the ceiling being exceeded.
pub(crate) fn checked_sum(values: &[usize], budget: NetworkBudget) -> Result<usize, NetworkError> {
    values
        .iter()
        .try_fold(0_usize, |sum, value| sum.checked_add(*value))
        .ok_or(NetworkError::MemoryLimit {
            required: usize::MAX,
            limit: budget.memory_limit_bytes,
        })
}

/// Refuses a workspace that would exceed the caller's ceiling.
pub(crate) fn check_memory(required: usize, budget: NetworkBudget) -> Result<(), NetworkError> {
    if required > budget.memory_limit_bytes {
        Err(NetworkError::MemoryLimit {
            required,
            limit: budget.memory_limit_bytes,
        })
    } else {
        Ok(())
    }
}
