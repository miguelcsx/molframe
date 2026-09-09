//! Owned results whose allocation charge follows their lifetime.

use super::{MemoryBudgetError, MemoryReservation};
use std::ops::Deref;

/// An owned result retaining the reservation for its backing allocations.
///
/// Sharing `Arc<Retained<T>>` shares one charge. The value is borrowed rather
/// than extracted, so an array view cannot silently outlive its accounting.
#[derive(Debug)]
pub struct Retained<T> {
    value: T,
    reservation: MemoryReservation,
}

impl<T> Retained<T> {
    /// Transfers an existing allocation reservation to its completed result.
    ///
    /// The caller must reserve before allocation and include every retained
    /// backing buffer in `bytes`. Unused scratch capacity is released here.
    ///
    /// # Errors
    ///
    /// Rejects an understated reservation without publishing the value.
    pub fn new(
        value: T,
        reservation: MemoryReservation,
        bytes: usize,
    ) -> Result<Self, MemoryBudgetError> {
        let mut retained = Self { value, reservation };
        if bytes > retained.reservation.bytes() {
            return Err(MemoryBudgetError::Exhausted {
                requested: bytes,
                available: retained.reservation.bytes(),
            });
        }
        retained.reservation.shrink_to(bytes);
        Ok(retained)
    }

    /// Charged bytes retained by this result.
    #[must_use]
    pub const fn retained_bytes(&self) -> usize {
        self.reservation.bytes()
    }
}

impl<T> Deref for Retained<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.value
    }
}

#[cfg(test)]
#[path = "retained_tests.rs"]
mod tests;
