//! The tables compiled into the library: none, so this layer is absent.

use crate::table_memory::MemoryTables;
use std::sync::Arc;

/// No layer: the library compiles in no table.
pub(crate) const fn bundled() -> Option<Arc<MemoryTables>> {
    None
}
