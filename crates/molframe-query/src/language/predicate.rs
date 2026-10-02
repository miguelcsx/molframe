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
    compare, is_residue_id_column, require_namespace, require_numeric, residue_symbol, same_by_key,
    symbol_residue_level, symbol_value,
};
use molframe_core::chunk::AtomChunk;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};
use scan::{SelectionCursor, scan_residues, selection_from_sorted, visit};
use std::collections::{BTreeSet, HashMap, hash_map::Entry};
use values::{numeric, residue_level, residue_text, text_resolved};

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
    let test = |mut actual: f64| {
        if absolute {
            actual = actual.abs();
        }
        compare(actual, expected, operator, tolerance)
    };
    // The two columns most often thresholded are read from each chunk directly.
    let chunked = match column {
        Column::BFactor => Some(scan::F32Column {
            plain: AtomChunk::b_factors_plain,
            validity: AtomChunk::b_factor_validity,
            get: AtomChunk::b_factor,
        }),
        Column::Occupancy => Some(scan::F32Column {
            plain: AtomChunk::occupancies_plain,
            validity: AtomChunk::occupancy_validity,
            get: AtomChunk::occupancy,
        }),
        _ => None,
    };
    if let Some(chunked) = chunked {
        return Ok(scan::scan_f32(structure, universe, &chunked, test));
    }
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
        return Ok(scan_residues(structure, universe, |_, residue| {
            residue_number(residue, column, policy)
                .is_some_and(|number| matcher.matches_parts(number, residue_insertion(residue)))
        }));
    }

    if column == Column::Element {
        return Ok(element_membership(structure, universe, values));
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

/// Selects atoms whose element is among `values`.
///
/// Elements are compared as elements, not as stored text: `Fe`, `FE` and `fe`
/// all name iron whatever case the file wrote, and a pattern such as `C*` is
/// matched against each element symbol without regard to case.
fn element_membership(
    structure: &Structure,
    universe: &AtomSelection,
    values: &[Box<str>],
) -> AtomSelection {
    let mut named = Vec::new();
    let mut patterns = Vec::new();
    for value in values {
        if value.bytes().any(|byte| matches!(byte, b'*' | b'?' | b'[')) {
            patterns.push(Glob::new(&value.to_ascii_uppercase()));
        } else if let Some(element) = molframe_core::Element::from_symbol(value) {
            named.push(element);
        }
    }
    if named.is_empty() && patterns.is_empty() {
        return AtomSelection::Empty;
    }
    // A named element is exactly what a chunk summary can rule out; a glob
    // pattern is not, so a pattern forces the full scan.
    let narrowed;
    let universe = if patterns.is_empty() {
        narrowed = universe.intersect(&super::prune::chunks_with_any_element(structure, &named));
        &narrowed
    } else {
        universe
    };
    scan(structure, universe, |context| {
        context.atom.element().is_some_and(|element| {
            named.contains(&element)
                || patterns
                    .iter()
                    .any(|pattern| pattern.matches(&element.symbol().to_ascii_uppercase()))
        })
    })
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
    if residue_level(column) {
        return scan_residues(structure, universe, |chain, residue| {
            let Some(actual) = residue_text(structure, chain, residue, column) else {
                return false;
            };
            *cache
                .entry(actual)
                .or_insert_with(|| globs.iter().any(|glob| glob.matches(actual)))
        });
    }
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
    symbols: &BTreeSet<molframe_core::symbol::SymbolId>,
) -> Result<AtomSelection, Diagnostic> {
    crate::annotation::require_available(structure, column)?;
    if symbols.is_empty() || universe.is_empty() {
        return Ok(AtomSelection::Empty);
    }
    if symbol_residue_level(column) {
        return Ok(scan_residues(structure, universe, |chain, residue| {
            residue_symbol(structure, chain, residue, column)
                .is_some_and(|value| symbols.contains(&value))
        }));
    }
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
        while let Some(position) = universe.next_before(atoms.end) {
            let Some(atom) = structure.atom(molframe_core::AtomIndex::new(position)) else {
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

#[cfg(test)]
#[path = "predicate_scan_tests.rs"]
mod scan_tests;
