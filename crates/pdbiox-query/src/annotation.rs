//! Access to conventional typed custom atom annotations.

use crate::ast::Column;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::structure::Structure;
use pdbiox_core::{AtomAnnotation, Presence, SymbolId};

pub(super) fn name(column: Column) -> Option<&'static str> {
    match column {
        Column::SegmentId => Some(pdbiox_core::SEGMENT_ID_ANNOTATION),
        Column::Charge => Some(pdbiox_core::PARTIAL_CHARGE_ANNOTATION),
        Column::Plddt => Some(pdbiox_core::PLDDT_ANNOTATION),
        Column::Pae => Some(pdbiox_core::PAE_ANNOTATION),
        _ => None,
    }
}

pub(super) fn require_available(structure: &Structure, column: Column) -> Result<(), Diagnostic> {
    let Some(name) = name(column) else {
        return Ok(());
    };
    let Some(annotation) = structure.annotations().get(name) else {
        return Err(Diagnostic::new(Code::E4003)
            .with_context("required", "requested annotation column")
            .with_context("annotation", name));
    };
    let compatible = if column == Column::SegmentId {
        matches!(annotation, AtomAnnotation::Symbol(_))
    } else {
        matches!(
            annotation,
            AtomAnnotation::Integer(_) | AtomAnnotation::Real(_)
        )
    };
    if compatible {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4002).with_context("annotation", name))
    }
}

pub(super) fn number(structure: &Structure, atom: u32, name: &str) -> Option<f64> {
    match structure.annotations().get(name)? {
        AtomAnnotation::Integer(column) => column
            .get(atom)
            .filter(|(_, presence)| *presence == Presence::Present)
            .map(|(value, _)| value as f64),
        AtomAnnotation::Real(column) => column
            .get(atom)
            .filter(|(_, presence)| *presence == Presence::Present)
            .map(|(value, _)| value),
        _ => None,
    }
}

pub(super) fn symbol(structure: &Structure, atom: u32, name: &str) -> Option<SymbolId> {
    let AtomAnnotation::Symbol(column) = structure.annotations().get(name)? else {
        return None;
    };
    column
        .get(atom)
        .filter(|(_, presence)| *presence == Presence::Present)
        .map(|(value, _)| value)
}
