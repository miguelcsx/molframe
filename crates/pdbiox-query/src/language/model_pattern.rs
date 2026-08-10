//! Model-level numeric membership over the structure universe.

use crate::ast::Column;
use crate::predicate_pattern::{NumericMatcher, NumericPattern};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;

/// Returns the whole universe when any model matches a compiled numeric pattern.
///
/// For `M` models and `P` patterns, runtime is `O(P log P + M log P)` for
/// multi-pattern inputs and `O(M)` for one pattern. The model index conversion
/// is checked because the query numeric domain must never receive a truncated
/// platform-sized index.
pub(super) fn model_membership(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    patterns: &[NumericPattern],
) -> AtomSelection {
    let matcher = NumericMatcher::from_patterns(patterns);
    let has_match = if column == Column::ModelIndex {
        (0..structure.model_count()).any(|model| {
            u32::try_from(model)
                .ok()
                .is_some_and(|model| matcher.matches(f64::from(model)))
        })
    } else {
        structure
            .data()
            .models()
            .filter_map(pdbiox_core::structure::ModelRef::number)
            .any(|number| matcher.matches(f64::from(number)))
    };

    if has_match {
        universe.clone()
    } else {
        AtomSelection::Empty
    }
}
