//! Column predicates, ranges, macros and same-column expansion.

use crate::ast::{Column, Operator};
use crate::glob::Glob;
use crate::predicate_pattern::{
    NumericMatcher, NumericPattern, ResidueMatcher, ResiduePattern, residue_insertion,
    residue_number,
};
use pdbiox_core::contract::{AnalysisPolicy, Namespace};
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};
use std::collections::{BTreeSet, HashMap};

#[path = "predicate_helpers.rs"]
mod helpers;
pub(crate) use helpers::atom_selector;
use helpers::{
    compare, entity_type, is_residue_id_column, require_namespace, require_numeric, same_by_key,
    symbol_value,
};

#[derive(Clone, Copy)]
pub(super) struct AtomContext<'a> {
    pub(super) atom: AtomRef<'a>,
    pub(super) residue: ResidueRef<'a>,
    pub(super) chain: ChainRef<'a>,
}

/// Selects atoms whose numeric `column` satisfies `operator expected`.
///
/// `absolute` applies `abs()` to the observed value before comparison. The
/// structure is scanned once, so runtime is `O(A * C_u)` where `A` is the
/// number of atoms and `C_u` is the cost of `universe.contains`; output space
/// is `O(M)` for `M` matching atom indices.
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
///
/// Numeric and residue pattern sets are compiled once, reducing repeated
/// membership checks to `O(log P)` for `P > 1`. Textual glob evaluation is
/// memoized per distinct borrowed value, avoiding per-key string allocation.
pub(crate) fn membership(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    values: &[Box<str>],
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_namespace(column, policy)?;
    crate::annotation::require_available(structure, column)?;

    if values.is_empty() {
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

/// Scans text-valued columns with a borrowed-value glob-result cache.
///
/// Each distinct textual value is evaluated against `globs` at most once.
/// Expected runtime is `O(A + U * G)` hash operations/glob evaluations for
/// `A` atoms, `U` distinct values, and `G` glob patterns, with `O(U)` cache
/// space and no per-key string allocation.
fn text_membership<'a>(
    structure: &'a Structure,
    universe: &AtomSelection,
    column: Column,
    globs: &[Glob],
    matches_empty_alternate: bool,
) -> AtomSelection {
    let mut cache = HashMap::<&'a str, bool>::new();
    let mut selected = Vec::new();

    visit(structure, |context| {
        let index = context.atom.index().get();

        if !universe.contains(index) {
            return;
        }

        let Some(actual) = text_resolved(structure, context, column) else {
            return;
        };

        if matches_empty_alternate && actual.is_empty() {
            selected.push(index);
            return;
        }

        let matched = if let Some(matched) = cache.get(actual) {
            *matched
        } else {
            let matched = globs.iter().any(|glob| glob.matches(actual));
            cache.insert(actual, matched);
            matched
        };

        if matched {
            selected.push(index);
        }
    });

    AtomSelection::from_sorted(selected)
}

/// Expands `selected` to every atom in `universe` sharing the same `column` value.
///
/// Numeric keys use their exact IEEE-754 bit representation, while textual
/// keys borrow directly from `structure`. Expected runtime is `O(A)` with
/// `O(K)` key space for `K` distinct selected values, without per-atom string
/// allocation.
pub(crate) fn same_column(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
    column: Column,
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_namespace(column, policy)?;
    crate::annotation::require_available(structure, column)?;

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
///
/// Element queries retain chunk-level pruning and store resolved elements in a
/// contiguous `Vec`. Runtime is proportional to visited candidate atoms plus
/// chunk pruning, with `O(M)` output space.
pub(crate) fn membership_symbols(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    symbols: &BTreeSet<pdbiox_core::symbol::SymbolId>,
) -> Result<AtomSelection, Diagnostic> {
    crate::annotation::require_available(structure, column)?;

    if symbols.is_empty() {
        return Ok(AtomSelection::Empty);
    }

    let elements: Vec<pdbiox_core::Element> = if column == Column::Element {
        symbols
            .iter()
            .filter_map(|symbol| structure.resolve(*symbol))
            .filter_map(pdbiox_core::Element::from_symbol)
            .collect()
    } else {
        Vec::new()
    };

    if column == Column::Element && elements.is_empty() {
        return Ok(AtomSelection::Empty);
    }

    let mut selected = Vec::new();

    for chunk in structure.data().chunks.iter() {
        if column == Column::Element
            && !elements
                .iter()
                .any(|element| chunk.stats().elements.contains(*element))
        {
            continue;
        }

        for position in chunk.atoms() {
            if !universe.contains(position) {
                continue;
            }

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

    Ok(AtomSelection::from_sorted(selected))
}

/// Scans atoms in structure order and returns those inside `universe` accepted by `accepts`.
///
/// Runtime is `O(A * C_u)` for `A` atoms and universe-membership cost `C_u`;
/// output space is `O(M)` for matching atoms.
#[inline]
pub(super) fn scan(
    structure: &Structure,
    universe: &AtomSelection,
    mut accepts: impl FnMut(AtomContext<'_>) -> bool,
) -> AtomSelection {
    let mut selected = Vec::new();

    visit(structure, |context| {
        let index = context.atom.index().get();

        if universe.contains(index) && accepts(context) {
            selected.push(index);
        }
    });

    AtomSelection::from_sorted(selected)
}

/// Visits every atom once while carrying its residue and chain context.
///
/// Runtime is `O(A)` and additional space is `O(1)` outside the caller's
/// closure state.
#[inline]
fn visit<'a>(structure: &'a Structure, mut visitor: impl FnMut(AtomContext<'a>)) {
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                visitor(AtomContext {
                    atom,
                    residue,
                    chain,
                });
            }
        }
    }
}

/// Extracts a numeric value for `column` from an atom context.
///
/// Returns `None` when the value or required policy mapping is unavailable.
/// Runtime is `O(1)` aside from annotation/property backend lookup costs and no
/// heap allocation is performed by this function.
#[inline]
fn numeric(
    structure: &Structure,
    context: AtomContext<'_>,
    column: Column,
    policy: &AnalysisPolicy,
) -> Option<f64> {
    let value = match column {
        Column::Index => f64::from(context.atom.index().get()),
        Column::ByNumber => f64::from(context.atom.index().get()) + 1.0,
        Column::AtomSiteId => f64::from(context.atom.atom_site_id()?),
        Column::ResidueIndex => f64::from(context.residue.index().get()),
        Column::ChainIndex => f64::from(context.chain.index().get()),
        Column::ModelIndex => 0.0,
        Column::Model => 1.0,
        Column::X => f64::from(context.atom.position()?[0]),
        Column::Y => f64::from(context.atom.position()?[1]),
        Column::Z => f64::from(context.atom.position()?[2]),
        Column::FormalCharge => match context.atom.formal_charge() {
            Some(charge) => f64::from(charge),
            None => crate::annotation::number(
                structure,
                context.atom.index().get(),
                pdbiox_core::FORMAL_CHARGE_ANNOTATION,
            )?,
        },
        Column::Mass => pdbiox_chem::element_properties(context.atom.element()?)?.atomic_weight,
        Column::Radius => f64::from(pdbiox_chem::vdw_radius(
            context.atom.element()?,
            radius_set(policy.vdw_radii)?,
        )?),
        Column::BFactor => f64::from(context.atom.b_factor()?),
        Column::Occupancy => f64::from(context.atom.occupancy()?),
        Column::Charge | Column::Plddt | Column::Pae => crate::annotation::number(
            structure,
            context.atom.index().get(),
            crate::annotation::name(column)?,
        )?,
        _ => return None,
    };

    Some(value)
}

/// Maps the analysis radii policy to the chemistry backend radius set.
///
/// Returns `None` for unsupported sets. Runtime and space are `O(1)`.
#[inline]
fn radius_set(set: pdbiox_core::contract::RadiiSet) -> Option<pdbiox_chem::RadiusSet> {
    Some(match set {
        pdbiox_core::contract::RadiiSet::Bondi => pdbiox_chem::RadiusSet::Bondi,
        pdbiox_core::contract::RadiiSet::AmberUnited => pdbiox_chem::RadiusSet::AmberUnited,
        pdbiox_core::contract::RadiiSet::Charmm => pdbiox_chem::RadiusSet::Charmm,
        pdbiox_core::contract::RadiiSet::Alvarez => pdbiox_chem::RadiusSet::Alvarez,
        _ => return None,
    })
}

/// Extracts a textual value after namespace resolution has already been performed.
///
/// The returned value borrows from `structure` or static storage, avoiding the
/// repeated namespace branch in hot scans.
#[inline]
fn text_resolved<'a>(
    structure: &'a Structure,
    context: AtomContext<'a>,
    column: Column,
) -> Option<&'a str> {
    match column {
        Column::LabelChain => context.chain.label(),
        Column::AuthChain => context.chain.auth_label(),
        Column::LabelResidueName => context.atom.component_name(),
        Column::AuthResidueName => context.residue.auth_name(),
        Column::LabelAtomName => context.atom.name(),
        Column::AuthAtomName => context.atom.auth_name(),
        Column::AlternateLocation => context
            .atom
            .alt_id()
            .and_then(pdbiox_core::symbol::AltId::symbol)
            .and_then(|symbol| structure.resolve(symbol))
            .or(Some("")),
        Column::Entity => context
            .chain
            .entity()
            .and_then(|entity| structure.data().topology.entities.id(entity))
            .and_then(|symbol| structure.resolve(symbol)),
        Column::EntityType => Some(entity_type(structure, context.chain)),
        Column::Element => context
            .atom
            .element()
            .map(pdbiox_core::element::Element::symbol),
        Column::InsertionCode => context.residue.ins_code().or(Some("")),
        Column::RecordType => Some(if context.residue.is_het() {
            "HETATM"
        } else {
            "ATOM"
        }),
        Column::SegmentId => crate::annotation::symbol(
            structure,
            context.atom.index().get(),
            pdbiox_core::SEGMENT_ID_ANNOTATION,
        )
        .and_then(|symbol| structure.resolve(symbol)),
        _ => None,
    }
}

/// Resolves policy-dependent logical columns to their concrete namespace columns.
///
/// Runtime and space are `O(1)`.
#[inline]
pub(crate) fn resolved_column(column: Column, policy: &AnalysisPolicy) -> Column {
    match (column, policy.identifiers) {
        (Column::Chain, Namespace::Label) => Column::LabelChain,
        (Column::Chain, Namespace::Auth) => Column::AuthChain,
        (Column::ResidueName, Namespace::Label) => Column::LabelResidueName,
        (Column::ResidueName, Namespace::Auth) => Column::AuthResidueName,
        (Column::AtomName, Namespace::Label) => Column::LabelAtomName,
        (Column::AtomName, Namespace::Auth) => Column::AuthAtomName,
        _ => column,
    }
}

#[cfg(test)]
#[path = "predicate_tests.rs"]
mod tests;
