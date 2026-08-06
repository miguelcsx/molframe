//! Writing the legacy fixed-column format, or refusing to.
//!
//! The format's fields have fixed widths, and a structure can outgrow every one
//! of them: a chain label of two characters, a serial past ninety-nine thousand,
//! a residue number past ten thousand, a coordinate too large for eight columns.
//!
//! Every one of those is checked before a single byte is written, and a
//! structure that does not fit is refused with a diagnostic naming the field and
//! offering a way out. Truncating a chain label from `AA` to `A` produces a file
//! that opens cleanly in any viewer and describes a different molecule, which is
//! the failure this library exists to prevent.

use crate::hybrid36;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::{Select, SelectAll};
use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};
use std::collections::HashMap;
use std::fmt::Write as _;

/// The largest magnitude an eight-column coordinate field can hold to three
/// decimal places.
const COORDINATE_LIMIT: f64 = 10_000.0;

/// How to write, and what to do about a structure that does not fit.
#[derive(Clone, Debug, Default)]
pub struct PdbOptions {
    /// Rename chains on the way out, so a structure with labels the format
    /// cannot hold can still be written deliberately.
    chain_map: HashMap<Box<str>, Box<str>>,
    /// Write serials and residue numbers past the decimal range in hybrid-36.
    ///
    /// Off by default: a consumer that does not implement the scheme would read
    /// a wrong number without noticing, so it is offered rather than assumed.
    hybrid36: bool,
}

impl PdbOptions {
    /// Options that refuse anything the format cannot hold.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Renames one chain on the way out.
    #[must_use]
    pub fn chain_map(mut self, from: impl Into<Box<str>>, to: impl Into<Box<str>>) -> Self {
        self.chain_map.insert(from.into(), to.into());
        self
    }

    /// Writes counts past the decimal range in hybrid-36 rather than refusing.
    #[must_use]
    pub const fn hybrid36(mut self, enabled: bool) -> Self {
        self.hybrid36 = enabled;
        self
    }

    /// The label a chain is written under.
    fn label_for<'a>(&'a self, label: &'a str) -> &'a str {
        match self.chain_map.get(label) {
            Some(replacement) => replacement,
            None => label,
        }
    }
}

/// Writes a structure, or explains why it cannot be written.
///
/// # Errors
///
/// Returns every capacity the structure exceeds, so a caller sees the whole
/// problem rather than fixing one field at a time.
pub fn write(structure: &Structure, options: &PdbOptions) -> Result<String, Vec<Diagnostic>> {
    write_selected(structure, options, &SelectAll)
}

/// Writes the parts of a structure a filter accepts.
///
/// # Errors
///
/// Returns every capacity the accepted parts exceed.
pub fn write_selected(
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
) -> Result<String, Vec<Diagnostic>> {
    let refusals = check_capacity(structure, options, select);
    if !refusals.is_empty() {
        return Err(refusals);
    }

    let data = structure.data();
    let mut out = String::with_capacity(structure.atom_count() as usize * 81);

    if let Some(id) = &data.entry.id {
        let _ = writeln!(out, "HEADER    {:<40}{:<9}   {:>4}", "", "", id);
    }
    if let Some(cell) = data.cell {
        let _ = writeln!(
            out,
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2}",
            cell.lengths[0],
            cell.lengths[1],
            cell.lengths[2],
            cell.angles[0],
            cell.angles[1],
            cell.angles[2],
        );
    }

    let mut serial = 1i64;
    for chain in data.chains() {
        if !select.accept_chain(chain.index()) {
            continue;
        }
        let label = chain_label(structure, &chain);
        let label = options.label_for(&label);

        for residue in chain.residues() {
            if !select.accept_residue(residue.index()) {
                continue;
            }
            for atom in residue.atoms() {
                if !select.accept_atom(atom.index()) {
                    continue;
                }
                write_atom(&mut out, structure, &atom, &residue, label, serial, options);
                serial += 1;
            }
        }
        let _ = writeln!(out, "TER   {:>5}", serial_field(serial, options));
        serial += 1;
    }
    out.push_str("END\n");
    Ok(out)
}

fn write_atom(
    out: &mut String,
    structure: &Structure,
    atom: &AtomRef<'_>,
    residue: &ResidueRef<'_>,
    chain: &str,
    serial: i64,
    options: &PdbOptions,
) {
    let record = if residue.is_het() { "HETATM" } else { "ATOM  " };
    let name = match atom.name() {
        Some(name) => name,
        None => "",
    };
    let element = match atom.element() {
        Some(element) => element.symbol(),
        None => "",
    };
    // An atom whose position was never recorded still has a row here, but the
    // format has no way to say "unrecorded" — so it is left out of the file
    // rather than written as an origin the experiment never observed.
    let Some(position) = atom.position() else {
        return;
    };
    let seq = match residue.auth_seq_id().or_else(|| residue.label_seq_id()) {
        Some(seq) => i64::from(seq),
        None => 0,
    };

    let _ = writeln!(
        out,
        "{record}{serial:>5} {name:<4}{alt:1}{comp:>3} {chain:>1}{seq}{ins:1}   \
         {x:8.3}{y:8.3}{z:8.3}{occ:6.2}{b:6.2}          {element:>2}",
        serial = serial_field(serial, options),
        name = atom_name_field(name),
        alt = alt_field(structure, atom),
        comp = match residue.name() {
            Some(comp) => comp,
            None => "",
        },
        seq = residue_field(seq, options),
        ins = match residue.ins_code() {
            Some(code) => code,
            None => " ",
        },
        x = f64::from(position[0]),
        y = f64::from(position[1]),
        z = f64::from(position[2]),
        occ = match atom.occupancy() {
            Some(occupancy) => f64::from(occupancy),
            None => 1.0,
        },
        b = match atom.b_factor() {
            Some(b_factor) => f64::from(b_factor),
            None => 0.0,
        },
    );
}

/// The four-column name field.
///
/// A name shorter than four characters starts in the second column, which is the
/// convention that lets a reader tell a one-letter element from a two-letter one.
fn atom_name_field(name: &str) -> String {
    if name.len() >= 4 {
        name.to_owned()
    } else {
        format!(" {name:<3}")
    }
}

/// The one-column alternate-location field.
///
/// A label longer than one character has already been refused by the capacity
/// check, so anything reaching here fits.
fn alt_field<'a>(structure: &'a Structure, atom: &AtomRef<'a>) -> &'a str {
    let Some(alt) = atom.alt_id() else { return " " };
    let Some(symbol) = alt.symbol() else {
        return " ";
    };
    match structure.resolve(symbol) {
        Some(label) => label,
        None => " ",
    }
}

fn serial_field(serial: i64, options: &PdbOptions) -> String {
    match options
        .hybrid36
        .then(|| hybrid36::encode(serial, 5))
        .flatten()
    {
        Some(field) => field,
        None => format!("{serial:>5}"),
    }
}

fn residue_field(seq: i64, options: &PdbOptions) -> String {
    match options.hybrid36.then(|| hybrid36::encode(seq, 4)).flatten() {
        Some(field) => field,
        None => format!("{seq:>4}"),
    }
}

fn chain_label(structure: &Structure, chain: &ChainRef<'_>) -> Box<str> {
    let symbol = chain.auth_asym_id().or_else(|| chain.label_asym_id());
    match symbol.and_then(|symbol| structure.resolve(symbol)) {
        Some(label) => label.into(),
        None => "".into(),
    }
}

/// Everything about the structure that the format cannot hold.
fn check_capacity(
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
) -> Vec<Diagnostic> {
    let data = structure.data();
    let mut refusals = Vec::new();

    for chain in data.chains() {
        if !select.accept_chain(chain.index()) {
            continue;
        }
        let original = chain_label(structure, &chain);
        let written = options.label_for(&original);
        if written.chars().count() > 1 {
            refusals.push(
                Diagnostic::new(Code::E4102)
                    .with_context("chain", original.to_string())
                    .with_context("written as", written.to_string()),
            );
        }
        for residue in chain.residues() {
            if !select.accept_residue(residue.index()) {
                continue;
            }
            check_residue(&residue, options, &mut refusals);
        }
    }

    let atoms = i64::from(structure.atom_count());
    if !options.hybrid36 && hybrid36::needs_hybrid36(atoms, 5) {
        refusals.push(
            Diagnostic::new(Code::E4101)
                .with_context("atoms", atoms.to_string())
                .with_context("field holds", "99999"),
        );
    }
    refusals
}

fn check_residue(residue: &ResidueRef<'_>, options: &PdbOptions, refusals: &mut Vec<Diagnostic>) {
    if let Some(seq) = residue.auth_seq_id().or_else(|| residue.label_seq_id())
        && !options.hybrid36
        && hybrid36::needs_hybrid36(i64::from(seq), 4)
    {
        refusals.push(
            Diagnostic::new(Code::E4103)
                .with_context("residue", seq.to_string())
                .with_context("field holds", "9999"),
        );
    }
    for atom in residue.atoms() {
        let Some(position) = atom.position() else {
            continue;
        };
        if position
            .iter()
            .any(|value| f64::from(*value).abs() >= COORDINATE_LIMIT)
        {
            refusals.push(
                Diagnostic::new(Code::E4104)
                    .with_context("atom", atom.index().to_string())
                    .with_context("position", format!("{position:?}")),
            );
        }
    }
}

#[cfg(test)]
#[path = "writer_tests.rs"]
mod tests;
