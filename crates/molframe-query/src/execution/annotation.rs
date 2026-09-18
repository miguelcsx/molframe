//! Access to conventional typed custom atom annotations.

use crate::ast::Column;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::Structure;
use molframe_core::{AtomAnnotation, Presence, SymbolId};

pub(crate) fn name(column: Column) -> Option<&'static str> {
    match column {
        Column::SegmentId => Some(molframe_core::SEGMENT_ID_ANNOTATION),
        Column::Charge => Some(molframe_core::PARTIAL_CHARGE_ANNOTATION),
        Column::Plddt => Some(molframe_core::PLDDT_ANNOTATION),
        Column::Pae => Some(molframe_core::PAE_ANNOTATION),
        _ => None,
    }
}

pub(crate) fn require_available(structure: &Structure, column: Column) -> Result<(), Diagnostic> {
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

pub(crate) fn number(structure: &Structure, atom: u32, name: &str) -> Option<f64> {
    match structure.annotations().get(name)? {
        AtomAnnotation::Integer(column) => column
            .get(atom)
            .filter(|(_, presence)| *presence == Presence::Present)
            .and_then(|(value, _)| exact_integer_as_f64(value)),
        AtomAnnotation::Real(column) => column
            .get(atom)
            .filter(|(_, presence)| *presence == Presence::Present)
            .map(|(value, _)| value),
        _ => None,
    }
}

/// Converts an integer annotation only when an `f64` preserves it exactly.
///
/// Numeric query literals use the binary64 domain. Treating a wider integer as
/// rounded would let a query select a neighboring value that was never stored.
#[inline]
fn exact_integer_as_f64(value: i64) -> Option<f64> {
    const MAX_EXACT_BINARY64_INTEGER: i64 = 1_i64 << f64::MANTISSA_DIGITS;
    const LOW_WORD_SCALE: f64 = 4_294_967_296.0;
    if !(-MAX_EXACT_BINARY64_INTEGER..=MAX_EXACT_BINARY64_INTEGER).contains(&value) {
        return None;
    }
    let high = i32::try_from(value.div_euclid(1_i64 << 32)).ok()?;
    let low = u32::try_from(value.rem_euclid(1_i64 << 32)).ok()?;
    Some(f64::from(high) * LOW_WORD_SCALE + f64::from(low))
}

pub(crate) fn symbol(structure: &Structure, atom: u32, name: &str) -> Option<SymbolId> {
    let AtomAnnotation::Symbol(column) = structure.annotations().get(name)? else {
        return None;
    };
    column
        .get(atom)
        .filter(|(_, presence)| *presence == Presence::Present)
        .map(|(value, _)| value)
}

pub(crate) fn boolean(structure: &Structure, atom: u32, name: &str) -> Option<bool> {
    let AtomAnnotation::Boolean(column) = structure.annotations().get(name)? else {
        return None;
    };
    column
        .get(atom)
        .filter(|(_, presence)| *presence == Presence::Present)
        .map(|(value, _)| value)
}

pub(crate) fn component_kind(
    structure: &Structure,
    atom: u32,
) -> Option<molframe_chem::ComponentKind> {
    let AtomAnnotation::Integer(column) = structure
        .annotations()
        .get(molframe_core::COMPONENT_KIND_ANNOTATION)?
    else {
        return None;
    };
    column
        .get(atom)
        .filter(|(_, presence)| *presence == Presence::Present)
        .and_then(|(code, _)| molframe_chem::ComponentKind::from_code(code))
}

pub(crate) fn polymer_atom_role(
    structure: &Structure,
    atom: u32,
) -> Option<molframe_chem::PolymerAtomRole> {
    let AtomAnnotation::Integer(column) = structure
        .annotations()
        .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION)?
    else {
        return None;
    };
    column
        .get(atom)
        .filter(|(_, presence)| *presence == Presence::Present)
        .and_then(|(code, _)| molframe_chem::PolymerAtomRole::from_code(code))
}

#[cfg(test)]
#[path = "annotation_tests.rs"]
mod tests;
