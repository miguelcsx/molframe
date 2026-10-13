//! Assembly of a decoded `StringArray` into a shared-dictionary column.
//!
//! One pass over the source strings builds the dictionary and one pass over the
//! row indices maps rows to it: `O(strings + rows)` expected time. An encoder
//! that already deduplicates its strings, which is the usual case, needs no
//! remapping table at all, so the only allocations are the dictionary, the
//! lookup and the row indices, plus one `Arc` per distinct string.

use super::DecodedStringColumn;
use super::decode::length_error;
use super::interner::Interner;
use molframe_core::diagnostic::{Code, Diagnostic};
use std::collections::HashMap;
use std::sync::Arc;

/// Builds a column from decoded row `indices` and string `offsets` into
/// `string_data`. A row index of `-1` is the empty string.
pub(super) fn build(
    indices: &[i64],
    offsets: &[i64],
    string_data: &str,
    mut interner: Option<&mut Interner>,
) -> Result<DecodedStringColumn, Diagnostic> {
    let source_count = offsets.len().saturating_sub(1);
    let mut lookup = HashMap::<&str, u32>::with_capacity(source_count);
    let mut dictionary = Vec::with_capacity(source_count);
    // Empty while every source string is distinct: its position is its index.
    let mut remap: Vec<u32> = Vec::new();
    let mut identity = true;
    for (position, pair) in offsets.windows(2).enumerate() {
        let (Ok(start), Ok(end)) = (usize::try_from(pair[0]), usize::try_from(pair[1])) else {
            return Err(length_error(offsets.len(), string_data.len()));
        };
        let Some(value) = string_data.get(start..end) else {
            return Err(length_error(end, string_data.len()));
        };
        let compact = intern(value, &mut lookup, &mut dictionary, interner.as_deref_mut())?;
        if identity && usize::try_from(compact) != Ok(position) {
            identity = false;
            remap.extend(0..compact_position(position)?);
        }
        if !identity {
            remap.push(compact);
        }
    }
    let mut output = Vec::with_capacity(indices.len());
    for (row, &index) in indices.iter().enumerate() {
        if index == -1 {
            output.push(intern(
                "",
                &mut lookup,
                &mut dictionary,
                interner.as_deref_mut(),
            )?);
            continue;
        }
        let Ok(index) = usize::try_from(index) else {
            return Err(length_error(row, source_count));
        };
        if index >= source_count {
            return Err(length_error(index, source_count));
        }
        output.push(if identity {
            compact_position(index)?
        } else {
            remap[index]
        });
    }
    Ok(DecodedStringColumn::from_validated_parts(
        dictionary, output,
    ))
}

fn compact_position(position: usize) -> Result<u32, Diagnostic> {
    u32::try_from(position).map_err(|_| too_many())
}

fn too_many() -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message("string dictionary exceeds u32 indices")
}

fn intern<'a>(
    value: &'a str,
    lookup: &mut HashMap<&'a str, u32>,
    dictionary: &mut Vec<Arc<str>>,
    interner: Option<&mut Interner>,
) -> Result<u32, Diagnostic> {
    if let Some(index) = lookup.get(value) {
        return Ok(*index);
    }
    let index = compact_position(dictionary.len())?;
    lookup.insert(value, index);
    dictionary.push(match interner {
        Some(interner) => interner.share(value),
        None => Arc::from(value),
    });
    Ok(index)
}
