//! Legacy PDB ranges resolve deposited endpoint identities in chain order.
//! Each record scans residues linearly with constant scratch space.

use crate::header::PdbHeaders;
use crate::{fixed, hybrid36};
use molframe_core::structure::{ResidueRef, StructureData};
use molframe_core::{SecondarySource, SecondaryStructure};

/// Lowers HELIX and SHEET records; every residue they name gets file provenance.
pub(super) fn read(
    data: &StructureData,
    headers: &PdbHeaders,
) -> (Vec<SecondaryStructure>, Vec<SecondarySource>) {
    let mut states = vec![SecondaryStructure::Unknown; data.topology.residues.len()];
    for (name, columns) in [
        ("HELIX", [20, 22, 25, 26, 32, 34, 37, 38]),
        ("SHEET", [22, 23, 26, 27, 33, 34, 37, 38]),
    ] {
        for record in headers.named(name) {
            let line = record.line();
            let kind = if name == "SHEET" {
                SecondaryStructure::Strand
            } else {
                helix_class(fixed::integer(line, 39, 40))
            };
            let begin = Endpoint::parse(line, &columns[..4]);
            let end = Endpoint::parse(line, &columns[4..]);
            if let (Some(begin), Some(end)) = (begin, end) {
                assign(data, &mut states, kind, begin, end);
            }
        }
    }
    let sources = states
        .iter()
        .map(|state| match state {
            SecondaryStructure::Unknown => SecondarySource::None,
            _ => SecondarySource::File,
        })
        .collect();
    (states, sources)
}

/// Classes 1, 3, 5 and 10 name alpha, pi, three-ten and polyproline helices.
pub(super) fn helix_class(class: Option<i64>) -> SecondaryStructure {
    match class {
        Some(1) => SecondaryStructure::AlphaHelix,
        Some(3) => SecondaryStructure::PiHelix,
        Some(5) => SecondaryStructure::ThreeTenHelix,
        Some(10) => SecondaryStructure::PolyProline,
        _ => SecondaryStructure::OtherHelix,
    }
}

#[derive(Clone, Copy)]
struct Endpoint<'a> {
    chain: &'a str,
    sequence: i32,
    insertion: &'a str,
}

impl<'a> Endpoint<'a> {
    fn parse(line: &'a str, columns: &[usize]) -> Option<Self> {
        Some(Self {
            chain: fixed::text(line, columns[0], columns[0]),
            sequence: i32::try_from(hybrid36::decode(
                fixed::raw(line, columns[1], columns[2]),
                4,
            )?)
            .ok()?,
            insertion: fixed::text(line, columns[3], columns[3]),
        })
    }

    fn matches(&self, residue: ResidueRef<'_>) -> bool {
        let insertion = match residue.ins_code() {
            Some(code) => code,
            None => "",
        };
        residue.auth_seq_id() == Some(self.sequence) && insertion == self.insertion
    }
}

fn assign(
    data: &StructureData,
    states: &mut [SecondaryStructure],
    kind: SecondaryStructure,
    begin: Endpoint<'_>,
    end: Endpoint<'_>,
) {
    if begin.chain != end.chain {
        return;
    }
    for chain in data.chains() {
        if chain.auth_label() != Some(begin.chain) {
            continue;
        }
        let mut first = None;
        let mut last = None;
        let mut ambiguous = false;
        for residue in chain.residues() {
            if begin.matches(residue) {
                ambiguous |= first.replace(residue.index().as_usize()).is_some();
            }
            if end.matches(residue) {
                ambiguous |= last.replace(residue.index().as_usize()).is_some();
            }
        }
        let (Some(first), Some(last)) = (first, last) else {
            continue;
        };
        if !ambiguous
            && first <= last
            && let Some(range) = states.get_mut(first..=last)
        {
            range.fill(kind);
        }
    }
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
