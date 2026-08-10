//! Per-atom fields carried by PDB-derived coordinate formats.

use super::ReadState;
use crate::reader::lines::Line;
use pdbiox_core::annotation::{
    ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn, AtomAnnotation,
    PARTIAL_CHARGE_ANNOTATION,
};
use pdbiox_core::column::Presence;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::Format;
use pdbiox_core::span::ByteSpan;
use pdbiox_core::symbol::SymbolId;

impl ReadState<'_> {
    pub(super) fn variant_atom(&mut self, line: &Line<'_>) {
        let mut fields = line.text.split_ascii_whitespace().rev();
        match self.variant {
            Format::Pqr => {
                let radius = fields.next().and_then(|value| value.parse::<f64>().ok());
                let charge = fields.next().and_then(|value| value.parse::<f64>().ok());
                self.partial_charges.push(real_annotation(charge));
                self.radii.push(real_annotation(radius));
                if charge.is_none() || radius.is_none() {
                    self.variant_error(line, "PQR charge or radius could not be read");
                }
            }
            Format::Pdbqt => {
                let atom_type = fields.next().filter(|value| !value.is_empty());
                let charge = fields.next().and_then(|value| value.parse::<f64>().ok());
                self.partial_charges.push(real_annotation(charge));
                match atom_type {
                    Some(atom_type) => {
                        let symbol = self.intern(atom_type);
                        self.autodock_types.push((symbol, Presence::Present));
                    }
                    None => self
                        .autodock_types
                        .push((SymbolId::from_raw(0), Presence::Unknown)),
                }
                if charge.is_none() || atom_type.is_none() {
                    self.variant_error(line, "PDBQT charge or atom type could not be read");
                }
            }
            _ => {}
        }
    }

    pub(super) fn finish_variant_annotations(&mut self) {
        if !self.partial_charges.is_empty() {
            let Ok(column) =
                AnnotationColumn::from_entries(std::mem::take(&mut self.partial_charges))
            else {
                self.findings.push(Diagnostic::new(Code::E6009));
                return;
            };
            let _ = self
                .data
                .annotations
                .insert(PARTIAL_CHARGE_ANNOTATION, AtomAnnotation::Real(column));
        }
        if !self.radii.is_empty() {
            let Ok(column) = AnnotationColumn::from_entries(std::mem::take(&mut self.radii)) else {
                self.findings.push(Diagnostic::new(Code::E6009));
                return;
            };
            let _ = self
                .data
                .annotations
                .insert(ATOM_RADIUS_ANNOTATION, AtomAnnotation::Real(column));
        }
        if !self.autodock_types.is_empty() {
            let Ok(column) =
                AnnotationColumn::from_entries(std::mem::take(&mut self.autodock_types))
            else {
                self.findings.push(Diagnostic::new(Code::E6009));
                return;
            };
            let _ = self
                .data
                .annotations
                .insert(AUTODOCK_TYPE_ANNOTATION, AtomAnnotation::Symbol(column));
        }
    }

    fn variant_error(&mut self, line: &Line<'_>, message: &'static str) {
        self.findings.push(
            Diagnostic::new(Code::E1202)
                .with_message(message)
                .at(ByteSpan::empty(line.at)),
        );
    }
}

fn real_annotation(value: Option<f64>) -> (f64, Presence) {
    match value {
        Some(value) => (value, Presence::Present),
        None => (0.0, Presence::Unknown),
    }
}
