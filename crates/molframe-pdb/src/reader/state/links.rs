//! `SSBOND` and `LINK` records as typed bonds between the atoms they name.
//!
//! These records state a deposited connection and its kind, so a structure read
//! from the legacy format keeps its disulfides and metal contacts in the bond
//! graph rather than only in the retained header text. Endpoints are named by
//! chain, residue number, insertion code and atom name; they are matched while
//! the atom rows stream past, and only for the identities some record asked for.

use super::ReadState;
use crate::fixed;
use crate::reader::lines::Line;
use molframe_core::bond::{BondOrder, BondProvenance, BondRecord};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::element::Element;
use molframe_core::index::AtomIndex;
use molframe_core::optional::OptionalSymbol;
use molframe_core::span::{ByteSpan, Position};
use std::collections::HashMap;
use std::fmt::Write as _;

/// Elements whose contacts a `LINK` record describes as coordination.
const METALS: [&str; 24] = [
    "LI", "NA", "K", "MG", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "MO", "CD", "HG", "PB", "SR",
    "BA", "AL", "CR", "V", "W", "PT", "AU", "AG",
];

/// One connection record awaiting its atoms.
pub(super) struct LinkRequest {
    first: String,
    second: String,
    disulfide: bool,
    at: Position,
}

/// What the stream resolved an identity to.
pub(super) type LinkAtoms = HashMap<String, Option<(AtomIndex, Element)>>;

fn identity(out: &mut String, name: &str, alt: &str, chain: &str, seq: &str, ins: &str) {
    out.clear();
    let _ = write!(out, "{name}|{alt}|{chain}|{seq}|{ins}");
}

impl ReadState<'_> {
    pub(super) fn ssbond(&mut self, line: &Line<'_>) {
        if !is_identity(line.text) {
            self.skip_image(line);
            return;
        }
        let first = Self::link_identity(line, "SG", "", (16, 16), (18, 21), (22, 22));
        let second = Self::link_identity(line, "SG", "", (30, 30), (32, 35), (36, 36));
        self.request(first, second, true, line);
    }

    pub(super) fn link(&mut self, line: &Line<'_>) {
        if !is_identity(line.text) {
            self.skip_image(line);
            return;
        }
        let first = Self::link_identity(
            line,
            fixed::text(line.text, 13, 16),
            fixed::text(line.text, 17, 17),
            (22, 22),
            (23, 26),
            (27, 27),
        );
        let second = Self::link_identity(
            line,
            fixed::text(line.text, 43, 46),
            fixed::text(line.text, 47, 47),
            (52, 52),
            (53, 56),
            (57, 57),
        );
        self.request(first, second, false, line);
    }

    fn skip_image(&mut self, line: &Line<'_>) {
        self.findings
            .push(Diagnostic::new(Code::W3013).at(ByteSpan::empty(line.at)));
    }

    fn link_identity(
        line: &Line<'_>,
        name: &str,
        alt: &str,
        chain: (usize, usize),
        seq: (usize, usize),
        ins: (usize, usize),
    ) -> String {
        let mut key = String::new();
        identity(
            &mut key,
            name,
            alt,
            fixed::text(line.text, chain.0, chain.1),
            fixed::text(line.text, seq.0, seq.1),
            fixed::text(line.text, ins.0, ins.1),
        );
        key
    }

    fn request(&mut self, first: String, second: String, disulfide: bool, line: &Line<'_>) {
        self.link_atoms.entry(first.clone()).or_insert(None);
        self.link_atoms.entry(second.clone()).or_insert(None);
        self.links.push(LinkRequest {
            first,
            second,
            disulfide,
            at: line.at,
        });
    }

    /// Records the atom just read when a connection record names its identity.
    pub(super) fn note_link_atom(&mut self, line: &Line<'_>, atom: AtomIndex, element: Element) {
        if self.link_atoms.is_empty() {
            return;
        }
        let mut key = std::mem::take(&mut self.key_scratch);
        let alt = fixed::text(line.text, 17, 17);
        identity(
            &mut key,
            fixed::text(line.text, 13, 16),
            alt,
            fixed::text(line.text, 22, 22),
            fixed::text(line.text, 23, 26),
            fixed::text(line.text, 27, 27),
        );
        if let Some(slot) = self.link_atoms.get_mut(&key)
            && slot.is_none()
        {
            *slot = Some((atom, element));
        }
        if !alt.is_empty() {
            identity(
                &mut key,
                fixed::text(line.text, 13, 16),
                "",
                fixed::text(line.text, 22, 22),
                fixed::text(line.text, 23, 26),
                fixed::text(line.text, 27, 27),
            );
            if let Some(slot) = self.link_atoms.get_mut(&key)
                && slot.is_none()
            {
                *slot = Some((atom, element));
            }
        }
        self.key_scratch = key;
    }

    /// The typed bonds every resolvable connection record stands for.
    pub(super) fn link_bonds(&mut self) -> Vec<BondRecord> {
        let requests = std::mem::take(&mut self.links);
        let mut bonds = Vec::with_capacity(requests.len());
        for request in requests {
            let first = self.link_atoms.get(&request.first).copied().flatten();
            let second = self.link_atoms.get(&request.second).copied().flatten();
            let (Some((atom_a, element_a)), Some((atom_b, element_b))) = (first, second) else {
                self.findings.push(
                    Diagnostic::new(Code::W3016)
                        .at(ByteSpan::empty(request.at))
                        .with_context("endpoint a", request.first)
                        .with_context("endpoint b", request.second),
                );
                continue;
            };
            let kind = if request.disulfide {
                "disulf"
            } else if is_metal(element_a) || is_metal(element_b) {
                "metalc"
            } else {
                "covale"
            };
            let kind = OptionalSymbol::some(self.intern(kind));
            bonds.push(BondRecord {
                atom_a,
                atom_b,
                order: BondOrder::Unknown,
                provenance: BondProvenance::File,
                kind,
            });
        }
        bonds
    }
}

fn is_metal(element: Element) -> bool {
    METALS.contains(&element.symbol().to_ascii_uppercase().as_str())
}

/// Both partners sit in the asymmetric unit when the symmetry operators read
/// as the identity, or are absent.
fn is_identity(line: &str) -> bool {
    [(60, 65), (67, 72)]
        .iter()
        .all(|&(from, to)| matches!(fixed::text(line, from, to), "" | "1555"))
}
