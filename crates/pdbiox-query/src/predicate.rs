//! Column predicates, ranges, macros and same-column expansion.

use crate::ast::{Column, Operator};
use crate::glob::Glob;
use crate::predicate_pattern::{NumericPattern, ResiduePattern, model_membership, residue_ordinal};
use pdbiox_core::contract::{AnalysisPolicy, Namespace};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};
use pdbiox_core::topology::EntityKind;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy)]
pub(super) struct AtomContext<'a> {
    pub(super) atom: AtomRef<'a>,
    pub(super) residue: ResidueRef<'a>,
    pub(super) chain: ChainRef<'a>,
}

pub(super) fn comparison(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    operator: Operator,
    expected: f64,
    absolute: bool,
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_numeric(structure, column)?;
    Ok(scan(structure, universe, |context| {
        numeric(structure, context, column, policy).is_some_and(|mut actual| {
            if absolute {
                actual = actual.abs();
            }
            compare(actual, expected, operator, policy)
        })
    }))
}

pub(super) fn membership(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    values: &[Box<str>],
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_namespace(column, policy)?;
    crate::annotation::require_available(structure, column)?;
    if column.is_numeric() {
        let patterns: Result<Vec<NumericPattern>, Diagnostic> = values
            .iter()
            .map(|value| NumericPattern::parse(value))
            .collect();
        let patterns = patterns?;
        if matches!(column, Column::Model | Column::ModelIndex) {
            return Ok(model_membership(structure, universe, column, &patterns));
        }
        return Ok(scan(structure, universe, |context| {
            numeric(structure, context, column, policy)
                .is_some_and(|actual| patterns.iter().any(|pattern| pattern.matches(actual)))
        }));
    }
    if matches!(
        column,
        Column::ResidueId | Column::LabelResidueId | Column::AuthResidueId
    ) {
        let patterns: Result<Vec<ResiduePattern>, Diagnostic> = values
            .iter()
            .map(|value| ResiduePattern::parse(value))
            .collect();
        let patterns = patterns?;
        return Ok(scan(structure, universe, |context| {
            residue_ordinal(context.residue, column, policy)
                .is_some_and(|actual| patterns.iter().any(|pattern| pattern.matches(&actual)))
        }));
    }

    let globs: Vec<Glob> = values.iter().map(|value| Glob::new(value)).collect();
    let mut cache: BTreeMap<Box<str>, bool> = BTreeMap::new();
    Ok(scan(structure, universe, |context| {
        let Some(actual) = text(structure, context, column, policy) else {
            return false;
        };
        if column == Column::AlternateLocation
            && values
                .iter()
                .any(|value| value.eq_ignore_ascii_case("none"))
            && actual.is_empty()
        {
            return true;
        }
        if let Some(matched) = cache.get(actual) {
            return *matched;
        }
        let matched = globs.iter().any(|glob| glob.matches(actual));
        cache.insert(actual.into(), matched);
        matched
    }))
}

pub(super) fn same_column(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
    column: Column,
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    require_namespace(column, policy)?;
    crate::annotation::require_available(structure, column)?;
    let mut keys = BTreeSet::new();
    visit(structure, |context| {
        if selected.contains(context.atom.index().get())
            && let Some(key) = key(structure, context, column, policy)
        {
            keys.insert(key);
        }
    });
    Ok(scan(structure, universe, |context| {
        key(structure, context, column, policy).is_some_and(|key| keys.contains(&key))
    }))
}

pub(super) fn membership_symbols(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    symbols: &BTreeSet<pdbiox_core::symbol::SymbolId>,
) -> Result<AtomSelection, Diagnostic> {
    crate::annotation::require_available(structure, column)?;
    let elements: BTreeSet<pdbiox_core::Element> = if column == Column::Element {
        symbols
            .iter()
            .filter_map(|symbol| structure.resolve(*symbol))
            .filter_map(pdbiox_core::Element::from_symbol)
            .collect()
    } else {
        BTreeSet::new()
    };
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

fn numeric(
    structure: &Structure,
    context: AtomContext<'_>,
    column: Column,
    policy: &AnalysisPolicy,
) -> Option<f64> {
    let value = match column {
        Column::Index => f64::from(context.atom.index().get()),
        Column::ByNumber => f64::from(context.atom.index().get().saturating_add(1)),
        Column::AtomSiteId => f64::from(context.atom.atom_site_id()?),
        Column::ResidueIndex => f64::from(context.residue.index().get()),
        Column::ChainIndex => f64::from(context.chain.index().get()),
        Column::ModelIndex => 0.0,
        Column::Model => 1.0,
        Column::X => f64::from(context.atom.position()?[0]),
        Column::Y => f64::from(context.atom.position()?[1]),
        Column::Z => f64::from(context.atom.position()?[2]),
        Column::FormalCharge => f64::from(context.atom.formal_charge()?),
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

fn radius_set(set: pdbiox_core::contract::RadiiSet) -> Option<pdbiox_chem::RadiusSet> {
    Some(match set {
        pdbiox_core::contract::RadiiSet::Bondi => pdbiox_chem::RadiusSet::Bondi,
        pdbiox_core::contract::RadiiSet::AmberUnited => pdbiox_chem::RadiusSet::AmberUnited,
        pdbiox_core::contract::RadiiSet::Charmm => pdbiox_chem::RadiusSet::Charmm,
        pdbiox_core::contract::RadiiSet::Alvarez => pdbiox_chem::RadiusSet::Alvarez,
        _ => return None,
    })
}

fn text<'a>(
    structure: &'a Structure,
    context: AtomContext<'a>,
    column: Column,
    policy: &AnalysisPolicy,
) -> Option<&'a str> {
    match resolved_column(column, policy) {
        Column::LabelChain => context.chain.label(),
        Column::AuthChain => context.chain.auth_label().or_else(|| context.chain.label()),
        Column::LabelResidueName => context.atom.component_name(),
        Column::AuthResidueName => context
            .residue
            .auth_name()
            .or_else(|| context.atom.component_name()),
        Column::LabelAtomName => context.atom.name(),
        Column::AuthAtomName => context.atom.auth_name().or_else(|| context.atom.name()),
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

fn symbol_value(
    structure: &Structure,
    context: AtomContext<'_>,
    column: Column,
) -> Option<pdbiox_core::symbol::SymbolId> {
    match column {
        Column::LabelChain => context.chain.label_asym_id(),
        Column::AuthChain => context
            .chain
            .auth_asym_id()
            .or_else(|| context.chain.label_asym_id()),
        Column::LabelResidueName => context.residue.label_comp_id(),
        Column::AuthResidueName => context
            .residue
            .auth_comp_id()
            .or_else(|| context.residue.label_comp_id()),
        Column::LabelAtomName => context.atom.name_symbol(),
        Column::AuthAtomName => context
            .atom
            .auth_name_symbol()
            .or_else(|| context.atom.name_symbol()),
        Column::AlternateLocation => context.atom.alt_id()?.symbol(),
        Column::Entity => context
            .chain
            .entity()
            .and_then(|entity| structure.data().topology.entities.id(entity)),
        Column::Element => structure
            .data()
            .dictionary
            .get(context.atom.element()?.symbol()),
        Column::InsertionCode => structure
            .data()
            .topology
            .residues
            .ins_code(context.residue.index()),
        Column::SegmentId => crate::annotation::symbol(
            structure,
            context.atom.index().get(),
            pdbiox_core::SEGMENT_ID_ANNOTATION,
        ),
        _ => None,
    }
}

fn require_namespace(column: Column, policy: &AnalysisPolicy) -> Result<(), Diagnostic> {
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

fn require_numeric(structure: &Structure, column: Column) -> Result<(), Diagnostic> {
    if column.is_numeric() {
        crate::annotation::require_available(structure, column)
    } else {
        Err(Diagnostic::new(Code::E4002).with_context("column", format!("{column:?}")))
    }
}

fn compare(actual: f64, expected: f64, operator: Operator, policy: &AnalysisPolicy) -> bool {
    let tolerance =
        policy.float_tolerance.absolute + policy.float_tolerance.relative * expected.abs();
    match operator {
        Operator::Less => actual < expected,
        Operator::LessEqual => actual <= expected,
        Operator::Greater => actual > expected,
        Operator::GreaterEqual => actual >= expected,
        Operator::Equal => (actual - expected).abs() <= tolerance,
        Operator::NotEqual => (actual - expected).abs() > tolerance,
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
enum ValueKey {
    Number(u64),
    Text(Box<str>),
}

fn key(
    structure: &Structure,
    context: AtomContext<'_>,
    column: Column,
    policy: &AnalysisPolicy,
) -> Option<ValueKey> {
    if column.is_numeric() {
        numeric(structure, context, column, policy).map(|value| ValueKey::Number(value.to_bits()))
    } else {
        text(structure, context, column, policy).map(|value| ValueKey::Text(value.into()))
    }
}

fn entity_type<'a>(structure: &'a Structure, chain: ChainRef<'_>) -> &'a str {
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

pub(super) fn atom_selector(
    structure: &Structure,
    universe: &AtomSelection,
    segment: &str,
    residue: i32,
    name: &str,
    policy: &AnalysisPolicy,
) -> Result<AtomSelection, Diagnostic> {
    crate::annotation::require_available(structure, Column::SegmentId)?;
    Ok(scan(structure, universe, |context| {
        let segment_matches = crate::annotation::symbol(
            structure,
            context.atom.index().get(),
            pdbiox_core::SEGMENT_ID_ANNOTATION,
        )
        .and_then(|symbol| structure.resolve(symbol))
            == Some(segment);
        let residue_matches = residue_ordinal(context.residue, Column::ResidueId, policy)
            .is_some_and(|ordinal| ordinal.number == residue && ordinal.insertion.is_empty());
        let name_matches = text(structure, context, Column::AtomName, policy) == Some(name);
        segment_matches && residue_matches && name_matches
    }))
}

#[cfg(test)]
#[path = "predicate_tests.rs"]
mod tests;
