//! Explicit legacy connectivity records.

use super::ReadState;
use crate::fixed;
use crate::hybrid36;
use crate::reader::lines::Line;
use molframe_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::AtomIndex;
use molframe_core::optional::OptionalSymbol;
use molframe_core::span::ByteSpan;
use std::collections::HashMap;

impl ReadState<'_> {
    pub(super) fn conect(&mut self, line: &Line<'_>) {
        let Some(source) = serial_field(line.text, 7, 11) else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("CONECT source serial could not be read")
                    .at(ByteSpan::empty(line.at)),
            );
            return;
        };
        for (from, to) in [(12, 16), (17, 21), (22, 26), (27, 31)] {
            if fixed::text(line.text, from, to).is_empty() {
                continue;
            }
            if let Some(target) = serial_field(line.text, from, to) {
                self.conect.push((source, target, line.at));
            } else {
                self.findings.push(
                    Diagnostic::new(Code::E1202)
                        .with_message("CONECT target serial could not be read")
                        .at(ByteSpan::empty(line.at)),
                );
            }
        }
    }

    pub(super) fn finish_bonds(&mut self) {
        if self.conect.is_empty() && self.links.is_empty() {
            return;
        }
        let mut records: Vec<BondRecord> = Vec::new();
        for (source, target, at) in &self.conect {
            let (Some(atom_a), Some(atom_b)) = (
                self.serial_to_atom.get(source),
                self.serial_to_atom.get(target),
            ) else {
                self.findings.push(
                    Diagnostic::new(Code::E3006)
                        .at(ByteSpan::empty(*at))
                        .with_context("serial a", source.to_string())
                        .with_context("serial b", target.to_string()),
                );
                continue;
            };
            records.push(BondRecord {
                atom_a: *atom_a,
                atom_b: *atom_b,
                order: BondOrder::Unknown,
                provenance: BondProvenance::File,
                // A CONECT record states that two atoms are bonded and nothing
                // about what kind of connection it is.
                kind: OptionalSymbol::NONE,
            });
        }
        let mut seen: HashMap<(AtomIndex, AtomIndex), usize> = records
            .iter()
            .enumerate()
            .map(|(position, bond)| (pair(bond), position))
            .collect();
        for bond in self.link_bonds() {
            // A connection that `CONECT` also names is one bond; the typed
            // record contributes only its kind.
            if let Some(&position) = seen.get(&pair(&bond)) {
                records[position].kind = bond.kind;
            } else {
                seen.insert(pair(&bond), records.len());
                records.push(bond);
            }
        }
        let mut bonds = BondTableBuilder::new();
        for bond in records {
            bonds.push(bond);
        }
        self.data.bonds = bonds.finish();
    }
}

fn pair(bond: &BondRecord) -> (AtomIndex, AtomIndex) {
    if bond.atom_a <= bond.atom_b {
        (bond.atom_a, bond.atom_b)
    } else {
        (bond.atom_b, bond.atom_a)
    }
}

fn serial_field(line: &str, from: usize, to: usize) -> Option<u32> {
    hybrid36::decode(fixed::raw(line, from, to), 5).and_then(|serial| u32::try_from(serial).ok())
}
