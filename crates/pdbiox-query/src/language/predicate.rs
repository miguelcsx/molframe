//! Column predicates, ranges, macros and same-column expansion.

#[path = "predicate_helpers.rs"]
mod helpers;
#[path = "predicate_scan.rs"]
mod scan;
#[path = "predicate_values.rs"]
mod values;

use crate::ast::{Column, Operator};
use crate::glob::Glob;
use crate::predicate_pattern::{
    NumericMatcher, NumericPattern, ResidueMatcher, ResiduePattern, residue_insertion,
    residue_number,
};
use helpers::{
    compare, is_residue_id_column, require_namespace, require_numeric, same_by_key, symbol_value,
};
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};
use scan::{SelectionCursor, selection_from_sorted, visit};
use std::collections::{BTreeSet, HashMap, hash_map::Entry};
use values::{numeric, text_resolved};

pub(crate) use helpers::atom_selector;
pub(crate) use scan::scan;
pub(crate) use values::resolved_column;

/// Atom-level query context with its owning residue and chain.
#[derive(Clone, Copy)]
pub(crate) struct AtomContext<'a> {
    pub(super) atom: AtomRef<'a>,
    pub(super) residue: ResidueRef<'a>,
    pub(super) chain: ChainRef<'a>,
}

/// Selects atoms whose numeric `column` satisfies `operator expected`.
pub(crate) fn comparison(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    operator: Operator,
    expected: f64,
    absolute: bool,
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_numeric(structure, column)?;
    let tolerance =
        policy.float_tolerance.absolute + policy.float_tolerance.relative * expected.abs();
    Ok(scan(structure, universe, |context| {
        numeric(structure, context, column, policy).is_some_and(|mut actual| {
            if absolute {
                actual = actual.abs();
            }
            compare(actual, expected, operator, tolerance)
        })
    }))
}

/// Selects atoms whose `column` matches any textual, numeric, or residue pattern.
pub(crate) fn membership(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    values: &[Box<str>],
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_namespace(column, policy)?;
    crate::annotation::require_available(structure, column)?;
    if values.is_empty() || universe.is_empty() {
        return Ok(AtomSelection::Empty);
    }

    if column.is_numeric() {
        let mut patterns = Vec::with_capacity(values.len());
        for value in values {
            patterns.push(NumericPattern::parse(value)?);
        }
        if matches!(column, Column::Model | Column::ModelIndex) {
            return Ok(crate::model_pattern::model_membership(
                structure, universe, column, &patterns,
            ));
        }
        let matcher = NumericMatcher::from_patterns(&patterns);
        return Ok(scan(structure, universe, |context| {
            numeric(structure, context, column, policy)
                .is_some_and(|actual| matcher.matches(actual))
        }));
    }

    if is_residue_id_column(column) {
        let mut patterns = Vec::with_capacity(values.len());
        for value in values {
            patterns.push(ResiduePattern::parse(value)?);
        }
        let matcher = ResidueMatcher::from_patterns(patterns);
        return Ok(scan(structure, universe, |context| {
            residue_number(context.residue, column, policy).is_some_and(|number| {
                matcher.matches_parts(number, residue_insertion(context.residue))
            })
        }));
    }

    let resolved = resolved_column(column, policy);
    let globs: Vec<Glob> = values.iter().map(|value| Glob::new(value)).collect();
    let matches_empty_alternate = column == Column::AlternateLocation
        && values
            .iter()
            .any(|value| value.eq_ignore_ascii_case("none"));
    Ok(text_membership(
        structure,
        universe,
        resolved,
        &globs,
        matches_empty_alternate,
    ))
}

fn text_membership<'a>(
    structure: &'a Structure,
    universe: &AtomSelection,
    column: Column,
    globs: &[Glob],
    matches_empty_alternate: bool,
) -> AtomSelection {
    if universe.is_empty() {
        return AtomSelection::Empty;
    }
    let mut cache = HashMap::<&'a str, bool>::new();
    let mut selected = Vec::new();
    visit(structure, universe, |context| {
        let Some(actual) = text_resolved(structure, context, column) else {
            return;
        };
        if matches_empty_alternate && actual.is_empty() {
            selected.push(context.atom.index().get());
            return;
        }
        let matched = match cache.entry(actual) {
            Entry::Occupied(entry) => *entry.get(),
            Entry::Vacant(entry) => {
                let matched = globs.iter().any(|glob| glob.matches(actual));
                entry.insert(matched);
                matched
            }
        };
        if matched {
            selected.push(context.atom.index().get());
        }
    });
    selection_from_sorted(selected)
}

/// Expands `selected` to every atom in `universe` sharing the same `column` value.
pub(crate) fn same_column(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
    column: Column,
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_namespace(column, policy)?;
    crate::annotation::require_available(structure, column)?;
    if universe.is_empty() || selected.is_empty() {
        return Ok(AtomSelection::Empty);
    }
    if column.is_numeric() {
        return Ok(same_by_key(structure, universe, selected, |context| {
            numeric(structure, context, column, policy).map(f64::to_bits)
        }));
    }
    let resolved = resolved_column(column, policy);
    Ok(same_by_key(structure, universe, selected, |context| {
        text_resolved(structure, context, resolved)
    }))
}

/// Selects atoms whose symbol-valued `column` is present in `symbols`.
pub(crate) fn membership_symbols(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    symbols: &BTreeSet<pdbiox_core::symbol::SymbolId>,
) -> Result<AtomSelection, Diagnostic> {
    crate::annotation::require_available(structure, column)?;
    if symbols.is_empty() || universe.is_empty() {
        return Ok(AtomSelection::Empty);
    }
    let elements = if column == Column::Element {
        let resolved: Vec<pdbiox_core::Element> = symbols
            .iter()
            .filter_map(|symbol| structure.resolve(*symbol))
            .filter_map(pdbiox_core::Element::from_symbol)
            .collect();
        if resolved.is_empty() {
            return Ok(AtomSelection::Empty);
        }
        Some(resolved)
    } else {
        None
    };

    let mut selected = Vec::new();
    let mut universe = SelectionCursor::new(universe.iter());
    for chunk in structure.data().chunks.iter() {
        if universe.is_exhausted() {
            break;
        }
        let atoms = chunk.atoms();
        universe.skip_before(atoms.start);
        if universe.is_exhausted() {
            break;
        }
        let chunk_matches = match &elements {
            Some(elements) => elements
                .iter()
                .any(|element| chunk.stats().elements.contains(*element)),
            None => true,
        };
        if !chunk_matches {
            universe.skip_before(atoms.end);
            continue;
        }
        while let Some(position) = universe.next_before(atoms.end) {
            let Some(atom) = structure.atom(pdbiox_core::AtomIndex::new(position)) else {
                continue;
            };
            let Some(residue) = atom.residue() else {
                continue;
            };
            let Some(chain_index) = structure
                .data()
                .topology
                .chains
                .containing(residue.index().get())
            else {
                continue;
            };
            let Some(chain) = structure.chain(chain_index) else {
                continue;
            };
            let context = AtomContext {
                atom,
                residue,
                chain,
            };
            if symbol_value(structure, context, column)
                .is_some_and(|value| symbols.contains(&value))
            {
                selected.push(position);
            }
        }
    }
    Ok(selection_from_sorted(selected))
}

#[cfg(test)]
#[path = "predicate_tests.rs"]
mod tests;
