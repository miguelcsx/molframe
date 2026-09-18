//! Numeric, textual, and namespace-resolved predicate values.

use super::AtomContext;
use crate::ast::Column;
use molframe_core::contract::{AnalysisPolicy, Namespace};
use molframe_core::structure::Structure;

pub(super) fn numeric(
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
                molframe_core::FORMAL_CHARGE_ANNOTATION,
            )?,
        },
        Column::Mass => molframe_chem::element_properties(context.atom.element()?)?.atomic_weight,
        Column::Radius => f64::from(molframe_chem::vdw_radius(
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

fn radius_set(set: molframe_core::contract::RadiiSet) -> Option<molframe_chem::RadiusSet> {
    Some(match set {
        molframe_core::contract::RadiiSet::Bondi => molframe_chem::RadiusSet::Bondi,
        molframe_core::contract::RadiiSet::AmberUnited => molframe_chem::RadiusSet::AmberUnited,
        molframe_core::contract::RadiiSet::Charmm => molframe_chem::RadiusSet::Charmm,
        molframe_core::contract::RadiiSet::Alvarez => molframe_chem::RadiusSet::Alvarez,
        _ => return None,
    })
}

pub(super) fn text_resolved<'a>(
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
            .and_then(molframe_core::symbol::AltId::symbol)
            .and_then(|symbol| structure.resolve(symbol))
            .or(Some("")),
        Column::Entity => context
            .chain
            .entity()
            .and_then(|entity| structure.data().topology.entities.id(entity))
            .and_then(|symbol| structure.resolve(symbol)),
        Column::EntityType => Some(super::helpers::entity_type(structure, context.chain)),
        Column::Element => context
            .atom
            .element()
            .map(molframe_core::element::Element::symbol),
        Column::InsertionCode => context.residue.ins_code().or(Some("")),
        Column::RecordType => Some(if context.residue.is_het() {
            "HETATM"
        } else {
            "ATOM"
        }),
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
