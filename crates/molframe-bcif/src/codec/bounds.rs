//! Declared-size ceilings that keep hostile codec chains from allocating.
//!
//! A chain states how many elements it expands to before a single byte is
//! decoded. The reader already knows how many rows the column may hold, so a
//! chain that claims more is rejected up front instead of being trusted with an
//! allocation.

use super::Encoding;
use molframe_core::diagnostic::{Code, Diagnostic};

/// Ceiling applied when a caller offers no row count of its own.
pub(super) const DEFAULT_ELEMENT_CEILING: usize = 1 << 28;

/// Rejects a chain whose declared sizes exceed `max_elements`.
///
/// # Errors
///
/// Returns `E1401` naming the offending declared size.
pub(crate) fn check_declared_sizes(
    chain: &[Encoding],
    max_elements: usize,
) -> Result<(), Diagnostic> {
    for encoding in chain {
        match encoding {
            Encoding::RunLength { src_size, .. } | Encoding::IntegerPacking { src_size, .. } => {
                if *src_size > max_elements {
                    return Err(Diagnostic::new(Code::E1401)
                        .with_message("declared element count exceeds the column bound")
                        .with_context("declared", src_size.to_string())
                        .with_context("bound", max_elements.to_string()));
                }
            }
            Encoding::StringArray {
                data_encoding,
                string_data,
                offset_encoding,
                ..
            } => {
                check_declared_sizes(data_encoding, max_elements)?;
                let offsets = max_elements
                    .saturating_add(string_data.len())
                    .saturating_add(1);
                check_declared_sizes(offset_encoding, offsets)?;
            }
            _ => {}
        }
    }
    Ok(())
}
