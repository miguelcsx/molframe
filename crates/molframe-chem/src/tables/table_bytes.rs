//! The encoded bytes of a table or a profile held in memory.
//!
//! A table that never touches a filesystem — compiled into the binary, fetched
//! by a browser, handed over by a host application — is still the same bytes a
//! cache file would hold, and decodes through the same [`TableKind`](crate::TableKind).
//! [`TableBytes`] is that payload: either borrowed from the binary for the
//! life of the process, or shared, so cloning it is a pointer copy and never a
//! copy of the table.

use core::fmt;
use std::sync::Arc;

/// Bytes of one table or profile, borrowed for the whole process or shared.
#[derive(Clone)]
pub struct TableBytes(Repr);

#[derive(Clone)]
enum Repr {
    Static(&'static [u8]),
    Shared(Arc<[u8]>),
}

impl TableBytes {
    /// Bytes that live as long as the process, such as an `include_bytes!` constant.
    #[must_use]
    pub const fn from_static(bytes: &'static [u8]) -> Self {
        Self(Repr::Static(bytes))
    }

    /// The bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        match &self.0 {
            Repr::Static(bytes) => bytes,
            Repr::Shared(bytes) => bytes,
        }
    }
}

impl AsRef<[u8]> for TableBytes {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl From<&'static [u8]> for TableBytes {
    fn from(bytes: &'static [u8]) -> Self {
        Self::from_static(bytes)
    }
}

impl From<&'static str> for TableBytes {
    fn from(text: &'static str) -> Self {
        Self::from_static(text.as_bytes())
    }
}

impl From<Arc<[u8]>> for TableBytes {
    fn from(bytes: Arc<[u8]>) -> Self {
        Self(Repr::Shared(bytes))
    }
}

impl From<Vec<u8>> for TableBytes {
    fn from(bytes: Vec<u8>) -> Self {
        Self(Repr::Shared(bytes.into()))
    }
}

impl From<String> for TableBytes {
    fn from(text: String) -> Self {
        Self::from(text.into_bytes())
    }
}

impl fmt::Debug for TableBytes {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        out.debug_struct("TableBytes")
            .field("len", &self.as_bytes().len())
            .finish()
    }
}
