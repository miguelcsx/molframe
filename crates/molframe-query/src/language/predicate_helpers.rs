//! Shared validation, comparison and expansion helpers for predicates.

use super::{AtomContext, resolved_column, scan, visit};
use crate::ast::{Column, Operator};
use crate::predicate_pattern::{residue_insertion, residue_number};
use molframe_core::contract::{AnalysisPolicy, Namespace};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{ChainRef, ResidueRef, Structure};
use molframe_core::topology::EntityKind;
use std::collections::HashSet;
use std::hash::Hash;

#[inline]
pub(super) fn require_namespace(column: Column, policy: &AnalysisPolicy) -> Result<(), Diagnostic> {
    if policy.identifiers == Namespace::Explicit
        && matches!(
            column,
            Column::Chain | Column::ResidueId | Column::ResidueName | Column::AtomName
        )
    {
        Err(Diagnostic::new(Code::E6001))
    } else {
        Ok(())
    }
}

pub(super) fn require_numeric(structure: &Structure, column: Column) -> Result<(), Diagnostic> {
    if column.is_numeric() {
        crate::annotation::require_available(structure, column)
    } else {
        Err(Diagnostic::new(Code::E4002).with_context("column", format!("{column:?}")))
    }
}

#[inline]
pub(super) fn compare(actual: f64, expected: f64, operator: Operator, tolerance: f64) -> bool {
    match operator {
        Operator::Less => actual < expected,
        Operator::LessEqual => actual <= expected,
        Operator::Greater => actual > expected,
        Operator::GreaterEqual => actual >= expected,
        Operator::Equal => (actual - expected).abs() <= tolerance,
        Operator::NotEqual => (actual - expected).abs() > tolerance,
    }
}

#[inline]
pub(super) fn is_residue_id_column(column: Column) -> bool {
    matches!(
        column,
        Column::ResidueId | Column::LabelResidueId | Column::AuthResidueId
    )
}

pub(super) fn same_by_key<'a, K>(
    structure: &'a Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
    mut key_for: impl FnMut(AtomContext<'a>) -> Option<K>,
) -> AtomSelection
where
    K: Eq + Hash,
{
    let mut keys = HashSet::new();
    visit(structure, selected, |context: AtomContext<'_>| {
        if let Some(key) = key_for(context) {
            keys.insert(key);
        }
    });
    if keys.is_empty() {
        return AtomSelection::Empty;
    }
    let mut matches = Vec::new();
    visit(structure, universe, |context: AtomContext<'_>| {
        let index = context.atom.index().get();
        if key_for(context).is_some_and(|key| keys.contains(&key)) {
            matches.push(index);
        }
    });
    AtomSelection::from_sorted(matches)
}

#[inline]
pub(super) fn entity_type<'a>(structure: &'a Structure, chain: ChainRef<'_>) -> &'a str {
    match chain
        .entity()
        .and_then(|entity| structure.data().topology.entities.kind(entity))
    {
        Some(EntityKind::Polymer) => "polymer",
        Some(EntityKind::NonPolymer) => "non-polymer",
        Some(EntityKind::Water) => "water",
        Some(EntityKind::Branched) => "branched",
        _ => "unknown",
    }
}

pub(crate) fn atom_selector(
    structure: &Structure,
    universe: &AtomSelection,
    segment: &str,
    residue: i32,
    name: &str,
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    crate::annotation::require_available(structure, Column::SegmentId)?;
    let Some(segment_symbol) = structure.data().dictionary.get(segment) else {
        return Ok(AtomSelection::Empty);
    };
    let Some(name_symbol) = structure.data().dictionary.get(name) else {
        return Ok(AtomSelection::Empty);
    };
    let atom_name_column = resolved_column(Column::AtomName, policy);
    Ok(scan(structure, universe, |context| {
        let segment_matches = crate::annotation::symbol(
            structure,
            context.atom.index().get(),
            molframe_core::SEGMENT_ID_ANNOTATION,
        ) == Some(segment_symbol);
        let residue_matches = residue_number(context.residue, Column::ResidueId, policy)
            .is_some_and(|number| {
                number == residue && residue_insertion(context.residue).is_empty()
            });
        let name_matches = symbol_value(structure, context, atom_name_column) == Some(name_symbol);
        segment_matches && residue_matches && name_matches
    }))
}

/// Whether the interned value of `column` is decided by a residue and its chain.
pub(super) const fn symbol_residue_level(column: Column) -> bool {
    matches!(
        column,
        Column::LabelChain
            | Column::AuthChain
            | Column::LabelResidueName
            | Column::AuthResidueName
            | Column::Entity
            | Column::InsertionCode
    )
}

/// The interned value of a residue-level `column`; `None` for any other column.
pub(super) fn residue_symbol(
    structure: &Structure,
    chain: ChainRef<'_>,
    residue: ResidueRef<'_>,
    column: Column,
) -> Option<molframe_core::symbol::SymbolId> {
    match column {
        Column::LabelChain => chain.label_asym_id(),
        Column::AuthChain => chain.auth_asym_id(),
        Column::LabelResidueName => residue.label_comp_id(),
        Column::AuthResidueName => residue.auth_comp_id(),
        Column::Entity => chain
            .entity()
            .and_then(|entity| structure.data().topology.entities.id(entity)),
        Column::InsertionCode => structure.data().topology.residues.ins_code(residue.index()),
        _ => None,
    }
}

#[inline]
pub(super) fn symbol_value(
    structure: &Structure,
    context: AtomContext<'_>,
    column: Column,
) -> Option<molframe_core::symbol::SymbolId> {
    if symbol_residue_level(column) {
        return residue_symbol(structure, context.chain, context.residue, column);
    }
    match column {
        Column::LabelAtomName => context.atom.name_symbol(),
        // Storage keeps an author atom name only where it differs from the
        // label name, so an atom without one is named by its label.
        Column::AuthAtomName => context
            .atom
            .auth_name_symbol()
            .or_else(|| context.atom.name_symbol()),
        Column::AlternateLocation => context.atom.alt_id()?.symbol(),
        Column::Element => structure
            .data()
            .dictionary
            .get(context.atom.element()?.symbol()),
        Column::SegmentId => crate::annotation::symbol(
            structure,
            context.atom.index().get(),
            molframe_core::SEGMENT_ID_ANNOTATION,
        ),
        _ => None,
    }
}
