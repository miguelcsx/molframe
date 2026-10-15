//! Numeric, textual, and namespace-resolved predicate values.

use super::AtomContext;
use crate::ast::Column;
use molframe_chem::{RadiusSource, TableRequest, TableStore};
use molframe_core::contract::{AnalysisPolicy, Namespace};
use molframe_core::diagnostic::Diagnostic;
use molframe_core::structure::{ChainRef, ResidueRef, Structure};

pub(super) fn numeric(
    structure: &Structure,
    context: AtomContext<'_>,
    column: Column,
    radii: Option<&RadiusSource>,
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
                molframe_core::FORMAL_CHARGE_ANNOTATION,
            )?,
        },
        Column::Mass => molframe_chem::element_properties(context.atom.element()?)?.atomic_weight,
        Column::Radius => f64::from(radii?.radius(context.atom.element()?)?),
        Column::BFactor => f64::from(context.atom.b_factor()?),
        Column::Occupancy => f64::from(context.atom.occupancy()?),
        Column::Charge | Column::Plddt | Column::Pae | Column::Assembly => {
            crate::annotation::number(
                structure,
                context.atom.index().get(),
                crate::annotation::name(column)?,
            )?
        }
        _ => return None,
    };
    Some(value)
}

/// The radius table a numeric `column` needs, resolved once per evaluation.
///
/// Only `radius` reads one: the table the policy's radius set names, looked up in
/// the cache. A missing table is the typed table error, never a silent default.
pub(super) fn radius_source(
    column: Column,
    policy: &AnalysisPolicy,
) -> Result<Option<RadiusSource>, Diagnostic> {
    if column != Column::Radius {
        return Ok(None);
    }
    let request = TableRequest::named(policy.vdw_radii.name());
    RadiusSource::resolve(&request, &TableStore::global())
        .map(Some)
        .map_err(Diagnostic::from)
}

/// Whether the text of `column` is decided by an atom's residue and chain alone.
///
/// Such a predicate needs one decision per residue rather than one per atom.
pub(super) const fn residue_level(column: Column) -> bool {
    matches!(
        column,
        Column::LabelChain
            | Column::AuthChain
            | Column::AuthResidueName
            | Column::Entity
            | Column::EntityType
            | Column::InsertionCode
            | Column::RecordType
    )
}

/// The text of a residue-level `column`; `None` for any other column.
pub(super) fn residue_text<'a>(
    structure: &'a Structure,
    chain: ChainRef<'a>,
    residue: ResidueRef<'a>,
    column: Column,
) -> Option<&'a str> {
    match column {
        Column::LabelChain => chain.label(),
        Column::AuthChain => chain.auth_label(),
        Column::AuthResidueName => residue.auth_name(),
        Column::Entity => chain
            .entity()
            .and_then(|entity| structure.data().topology.entities.id(entity))
            .and_then(|symbol| structure.resolve(symbol)),
        Column::EntityType => Some(super::columns::entity_type(structure, chain)),
        Column::InsertionCode => residue.ins_code().or(Some("")),
        Column::RecordType => Some(if residue.is_het() { "HETATM" } else { "ATOM" }),
        _ => None,
    }
}

pub(super) fn text_resolved<'a>(
    structure: &'a Structure,
    context: AtomContext<'a>,
    column: Column,
) -> Option<&'a str> {
    if residue_level(column) {
        return residue_text(structure, context.chain, context.residue, column);
    }
    match column {
        Column::LabelResidueName => context.atom.component_name(),
        Column::LabelAtomName => context.atom.name(),
        // Storage keeps an author atom name only where it differs from the
        // label name, so an atom without one is named by its label.
        Column::AuthAtomName => context.atom.auth_name().or_else(|| context.atom.name()),
        Column::AlternateLocation => context
            .atom
            .alt_id()
            .and_then(molframe_core::symbol::AltId::symbol)
            .and_then(|symbol| structure.resolve(symbol))
            .or(Some("")),
        Column::Element => context
            .atom
            .element()
            .map(molframe_core::element::Element::symbol),
        Column::SegmentId => crate::annotation::symbol(
            structure,
            context.atom.index().get(),
            molframe_core::SEGMENT_ID_ANNOTATION,
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
