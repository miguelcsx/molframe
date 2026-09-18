//! Per-atom fields carried by PDB-derived coordinate formats.

use super::ReadState;
use crate::reader::lines::Line;
use molframe_core::annotation::{
    ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn, AtomAnnotation,
    PARTIAL_CHARGE_ANNOTATION,
};
use molframe_core::column::Presence;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::Format;
use molframe_core::span::ByteSpan;
use molframe_core::symbol::SymbolId;

impl ReadState<'_> {
    pub(super) fn variant_atom(&mut self, line: &Line<'_>) {
        let fields = variant_fields(line, self.variant);
        match self.variant {
            Format::Pqr => {
                self.partial_charges.push(real_annotation(fields.charge));
                self.radii.push(real_annotation(fields.radius));
                if fields.charge.is_none() || fields.radius.is_none() {
                    self.variant_error(line, "PQR charge or radius could not be read");
                }
            }
            Format::Pdbqt => {
                self.partial_charges.push(real_annotation(fields.charge));
                match fields.atom_type {
                    Some(atom_type) => {
                        let symbol = self.intern(atom_type);
                        self.autodock_types.push((symbol, Presence::Present));
                    }
                    None => self
                        .autodock_types
                        .push((SymbolId::from_raw(0), Presence::Unknown)),
                }
                if fields.charge.is_none() || fields.atom_type.is_none() {
                    self.variant_error(line, "PDBQT charge or atom type could not be read");
                }
            }
            _ => {}
        }
    }

    pub(super) fn variant_atom_matches(&mut self, line: &Line<'_>) -> bool {
        let fields = variant_fields(line, self.variant);
        let position = self.atom_position as usize;
        match self.variant {
            Format::Pqr => {
                if fields.charge.is_none() || fields.radius.is_none() {
                    self.variant_error(line, "PQR charge or radius could not be read");
                }
                self.partial_charges.get(position) == Some(&real_annotation(fields.charge))
                    && self.radii.get(position) == Some(&real_annotation(fields.radius))
            }
            Format::Pdbqt => {
                if fields.charge.is_none() || fields.atom_type.is_none() {
                    self.variant_error(line, "PDBQT charge or atom type could not be read");
                }
                let charge_matches =
                    self.partial_charges.get(position) == Some(&real_annotation(fields.charge));
                let type_matches = match self.autodock_types.get(position) {
                    Some((symbol, Presence::Present)) => {
                        fields.atom_type.is_some_and(|atom_type| {
                            self.data.dictionary.resolve(*symbol) == Some(atom_type)
                        })
                    }
                    Some((_, presence)) => fields.atom_type.is_none() && !presence.is_present(),
                    None => false,
                };
                charge_matches && type_matches
            }
            _ => true,
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

struct VariantFields<'a> {
    charge: Option<f64>,
    radius: Option<f64>,
    atom_type: Option<&'a str>,
}

fn variant_fields<'a>(line: &'a Line<'_>, variant: Format) -> VariantFields<'a> {
    let mut fields = line.text.split_ascii_whitespace().rev();
    match variant {
        Format::Pqr => VariantFields {
            radius: fields.next().and_then(|value| value.parse::<f64>().ok()),
            charge: fields.next().and_then(|value| value.parse::<f64>().ok()),
            atom_type: None,
        },
        Format::Pdbqt => VariantFields {
            atom_type: fields.next().filter(|value| !value.is_empty()),
            charge: fields.next().and_then(|value| value.parse::<f64>().ok()),
            radius: None,
        },
        _ => VariantFields {
            charge: None,
            radius: None,
            atom_type: None,
        },
    }
}

fn real_annotation(value: Option<f64>) -> (f64, Presence) {
    match value {
        Some(value) => (value, Presence::Present),
        None => (0.0, Presence::Unknown),
    }
}
