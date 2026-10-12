//! Connectivity records: `SSBOND`, `LINK`, `CONECT`, and the `SCALE` matrix.
//!
//! A structure read from mmCIF carries its connections as bonds with a deposited
//! kind, and its cell as lengths and angles. Writing it to the legacy format
//! without these records leaves a file that opens cleanly and has lost every
//! disulfide, metal coordination and covalent link, so they are generated here
//! unless the structure's own raw headers already carry them.

use super::pdb::{
    PdbOptions, alt_field, atom_name_field, insertion_code, residue_field, serial_field,
};
use super::pdb::{atom_name, chain_label, component_name, residue_sequence};
use molframe_core::bond::BondProvenance;
use molframe_core::index::AtomIndex;
use molframe_core::io::Select;
use molframe_core::structure::{AtomRef, Structure, UnitCell};
use std::collections::{BTreeMap, HashMap};
use std::fmt;

/// Which legacy record a bond kind is written as.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Class {
    Disulfide,
    Link,
}

fn classify(structure: &Structure, bond: &molframe_core::bond::BondRecord) -> Option<Class> {
    if bond.provenance == BondProvenance::InferredDistance {
        return None;
    }
    let symbol = bond.kind.get()?;
    match structure.resolve(symbol)? {
        "disulf" => Some(Class::Disulfide),
        "covale" | "metalc" | "modres" => Some(Class::Link),
        text if text.starts_with("covale") => Some(Class::Link),
        _ => None,
    }
}

/// True when some bond is a bare `CONECT` statement.
pub(super) fn has_connect(structure: &Structure) -> bool {
    structure.data().bonds.iter().any(is_connect)
}

fn is_connect(bond: molframe_core::bond::BondRecord) -> bool {
    bond.provenance == BondProvenance::File && bond.kind.get().is_none()
}

struct Endpoint<'a> {
    name: &'a str,
    alt: &'a str,
    component: &'a str,
    chain: String,
    seq: String,
    ins: &'a str,
    position: Option<[f32; 3]>,
}

fn endpoint<'a>(
    structure: &'a Structure,
    atom: AtomRef<'a>,
    options: &'a PdbOptions,
) -> Option<Endpoint<'a>> {
    let namespace = options.identifier_namespace();
    let residue = atom.residue()?;
    let chain = structure
        .data()
        .topology
        .chains
        .containing(residue.index().get())?;
    let chain = structure.chain(chain)?;
    let label = options.label_for(chain_label(&chain, namespace)).to_owned();
    let seq = residue_sequence(residue, namespace)?;
    Some(Endpoint {
        name: atom_name(atom, namespace)?,
        alt: alt_field(structure, &atom),
        component: component_name(atom, residue, namespace)?,
        chain: label,
        seq: residue_field(i64::from(seq), options),
        ins: insertion_code(residue),
        position: atom.position(),
    })
}

fn distance(a: &Endpoint<'_>, b: &Endpoint<'_>) -> Option<f64> {
    let (a, b) = (a.position?, b.position?);
    let square: f64 = (0..3).map(|i| f64::from(a[i] - b[i]).powi(2)).sum();
    Some(square.sqrt())
}

/// Writes `SSBOND` and `LINK` for the deposited connections that have no raw
/// record in the structure's headers.
pub(super) fn write_links(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
    has_ssbond: bool,
    has_link: bool,
) {
    let mut ss = 0u32;
    for bond in structure.data().bonds.iter() {
        let Some(class) = classify(structure, &bond) else {
            continue;
        };
        if (class == Class::Disulfide && has_ssbond) || (class == Class::Link && has_link) {
            continue;
        }
        let (Some(a), Some(b)) = (structure.atom(bond.atom_a), structure.atom(bond.atom_b)) else {
            continue;
        };
        if !select.accept_atom(a.index()) || !select.accept_atom(b.index()) {
            continue;
        }
        let (Some(first), Some(second)) = (
            endpoint(structure, a, options),
            endpoint(structure, b, options),
        ) else {
            continue;
        };
        let length = distance(&first, &second).unwrap_or(0.0);
        match class {
            Class::Disulfide => {
                ss += 1;
                let _ = writeln!(
                    out,
                    "SSBOND {ss:>3} {:>3} {:1} {}{:1}   {:>3} {:1} {}{:1}                         \
                     1555   1555{length:>6.2}",
                    first.component,
                    first.chain,
                    first.seq,
                    first.ins,
                    second.component,
                    second.chain,
                    second.seq,
                    second.ins,
                );
            }
            Class::Link => {
                let _ = writeln!(
                    out,
                    "LINK        {}{:1}{:>3} {:1}{}{:1}               {}{:1}{:>3} {:1}{}{:1}  \
                     {:>6} {:>6} {length:5.2}",
                    atom_name_field(first.name),
                    first.alt,
                    first.component,
                    first.chain,
                    first.seq,
                    first.ins,
                    atom_name_field(second.name),
                    second.alt,
                    second.component,
                    second.chain,
                    second.seq,
                    second.ins,
                    "1555",
                    "1555",
                );
            }
        }
    }
}

/// Writes `CONECT` for the bare connectivity statements, both directions, four
/// partners to a line.
pub(super) fn write_connect(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
    serials: &HashMap<u32, i64>,
) {
    let mut partners: BTreeMap<i64, Vec<i64>> = BTreeMap::new();
    for bond in structure.data().bonds.iter().filter(|b| is_connect(*b)) {
        if !select.accept_atom(AtomIndex::new(bond.atom_a.get()))
            || !select.accept_atom(AtomIndex::new(bond.atom_b.get()))
        {
            continue;
        }
        let (Some(a), Some(b)) = (
            serials.get(&bond.atom_a.get()),
            serials.get(&bond.atom_b.get()),
        ) else {
            continue;
        };
        partners.entry(*a).or_default().push(*b);
        partners.entry(*b).or_default().push(*a);
    }
    for (source, mut targets) in partners {
        targets.sort_unstable();
        targets.dedup();
        for chunk in targets.chunks(4) {
            let _ = write!(out, "CONECT{}", serial_field(source, options));
            for target in chunk {
                let _ = write!(out, "{}", serial_field(*target, options));
            }
            let _ = writeln!(out);
        }
    }
}

/// Writes the three `SCALE` records, the inverse of the orthogonalization
/// matrix, so a reader can fractionalize without recomputing it.
pub(super) fn write_scale(out: &mut impl fmt::Write, cell: &UnitCell) {
    let [a, b, c] = cell.lengths;
    let [alpha, beta, gamma] = cell.angles.map(f64::to_radians);
    let sin_gamma = gamma.sin();
    if a <= 0.0 || b <= 0.0 || c <= 0.0 || sin_gamma.abs() < 1e-9 {
        return;
    }
    let volume_factor = (1.0 - alpha.cos().powi(2) - beta.cos().powi(2) - gamma.cos().powi(2)
        + 2.0 * alpha.cos() * beta.cos() * gamma.cos())
    .max(0.0)
    .sqrt();
    if volume_factor < 1e-9 {
        return;
    }
    let (o11, o12, o13) = (a, b * gamma.cos(), c * beta.cos());
    let o22 = b * sin_gamma;
    let o23 = c * (alpha.cos() - beta.cos() * gamma.cos()) / sin_gamma;
    let o33 = c * volume_factor / sin_gamma;
    let rows = [
        [
            1.0 / o11,
            -o12 / (o11 * o22),
            (o12 * o23 - o13 * o22) / (o11 * o22 * o33),
        ],
        [0.0, 1.0 / o22, -o23 / (o22 * o33)],
        [0.0, 0.0, 1.0 / o33],
    ];
    for (index, row) in rows.iter().enumerate() {
        let _ = writeln!(
            out,
            "SCALE{}    {:10.6}{:10.6}{:10.6}     {:10.5}",
            index + 1,
            row[0],
            row[1],
            row[2],
            0.0,
        );
    }
}
