//! Writing, under two contracts that are stated rather than implied.
//!
//! **Canonical** produces a valid, self-consistent file. It does not promise to
//! reproduce the original's formatting, and it does not carry categories the
//! structure does not model.
//!
//! **Preserving** writes a document back out with its categories, its items and
//! their order intact, and with every untouched value written the way it was
//! found — including whether it was quoted.
//!
//! Neither promises a byte-exact round trip after the data has been changed.
//! That promise cannot be kept honestly once a value differs in width, and a
//! library that implies it will be trusted where it should not be.

use crate::document::{CifValue, Document};
use crate::lexer::Quoting;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use std::fmt::Write as _;

/// Writes a document back out, preserving what it held.
///
/// Categories, items and their order survive, and so does the quoting of every
/// value, so a category this library has no interpretation for comes back
/// unharmed.
#[must_use]
pub fn write_preserving(document: &Document) -> String {
    let mut out = String::new();
    for block in document.blocks() {
        let _ = writeln!(out, "data_{}", block.name());
        for category in block.categories() {
            write_category(&mut out, category);
        }
    }
    out
}

fn write_category(out: &mut String, category: &crate::document::Category) {
    let rows = category.row_count();
    if rows == 1 {
        for item in category.items() {
            let Some(value) = category.value(item, 0) else {
                continue;
            };
            let quoting = category.column(item).and_then(|column| column.quoting(0));
            let _ = writeln!(
                out,
                "_{}.{:<30} {}",
                category.name(),
                item,
                render(value, quoting)
            );
        }
        out.push_str("#\n");
        return;
    }

    out.push_str("loop_\n");
    for item in category.items() {
        let _ = writeln!(out, "_{}.{item}", category.name());
    }
    for row in 0..rows {
        let mut first = true;
        for item in category.items() {
            let Some(value) = category.value(item, row) else {
                continue;
            };
            let quoting = category.column(item).and_then(|column| column.quoting(row));
            if !first {
                out.push(' ');
            }
            first = false;
            out.push_str(&render(value, quoting));
        }
        out.push('\n');
    }
    out.push_str("#\n");
}

/// Renders one value the way it was written, quoting it if it needs quoting.
fn render(value: &CifValue, quoting: Option<Quoting>) -> String {
    let text = match value {
        CifValue::Inapplicable => return ".".to_owned(),
        CifValue::Unknown => return "?".to_owned(),
        CifValue::Integer(number) => return number.to_string(),
        CifValue::Float(number) => return format!("{number}"),
        CifValue::Text(text) => text,
    };
    match quoting {
        Some(Quoting::Text) => format!("\n;{text}\n;"),
        Some(Quoting::Double) => format!("\"{text}\""),
        Some(Quoting::Single) => format!("'{text}'"),
        _ if needs_quoting(text) => format!("'{text}'"),
        _ => text.to_string(),
    }
}

/// Whether a bare value would be read back as something other than itself.
fn needs_quoting(text: &str) -> bool {
    text.is_empty()
        || text.chars().any(char::is_whitespace)
        || matches!(text, "." | "?")
        || text.starts_with(['_', '#', '\'', '"', '[', ']', '$'])
        || text.eq_ignore_ascii_case("loop_")
}

/// Writes a structure as a valid, self-consistent file.
///
/// This is a lowering in reverse: it carries what the structure models, which is
/// meaning rather than formatting.
#[must_use]
pub fn write_canonical(structure: &Structure) -> String {
    let data = structure.data();
    let name = match &data.entry.id {
        Some(id) => id.to_string(),
        None => "pdbiox".to_owned(),
    };
    let mut out = String::with_capacity(structure.atom_count() as usize * 100);
    let _ = writeln!(out, "data_{name}");
    out.push_str("#\n");

    if let Some(id) = &data.entry.id {
        let _ = writeln!(out, "_entry.id   {id}\n#");
    }
    if let Some(cell) = data.cell {
        let _ = writeln!(
            out,
            "_cell.length_a      {:.4}\n_cell.length_b      {:.4}\n_cell.length_c      {:.4}\n\
             _cell.angle_alpha   {:.4}\n_cell.angle_beta    {:.4}\n_cell.angle_gamma   {:.4}\n#",
            cell.lengths[0],
            cell.lengths[1],
            cell.lengths[2],
            cell.angles[0],
            cell.angles[1],
            cell.angles[2],
        );
    }

    out.push_str(ATOM_SITE_HEADER);
    let mut serial = 1u32;
    for (position, chain) in data.chains().enumerate() {
        let chain_label = normalised_label(position);
        let auth_label = label(
            structure,
            chain.auth_asym_id().or_else(|| chain.label_asym_id()),
        );
        for residue in chain.residues() {
            for atom in residue.atoms() {
                write_atom(
                    &mut out,
                    structure,
                    &atom,
                    &residue,
                    &chain_label,
                    &auth_label,
                    serial,
                );
                serial += 1;
            }
        }
    }
    out.push_str("#\n");
    out
}

/// The item names of the coordinate category, in the order they are written.
const ATOM_SITE_HEADER: &str = "loop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_alt_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
_atom_site.occupancy\n\
_atom_site.B_iso_or_equiv\n\
_atom_site.auth_seq_id\n\
_atom_site.auth_asym_id\n\
_atom_site.pdbx_PDB_model_num\n";

fn write_atom(
    out: &mut String,
    structure: &Structure,
    atom: &AtomRef<'_>,
    residue: &ResidueRef<'_>,
    chain: &str,
    auth_chain: &str,
    serial: u32,
) {
    let group = if residue.is_het() { "HETATM" } else { "ATOM" };
    let element = match atom.element() {
        Some(element) => element.symbol().to_uppercase(),
        None => "?".to_owned(),
    };
    let name = match atom.name() {
        Some(name) => name,
        None => "?",
    };
    let comp = match residue.name() {
        Some(comp) => comp,
        None => "?",
    };
    // An atom whose position was never recorded is written as unrecorded rather
    // than as an origin the experiment never observed.
    let [x, y, z] = match atom.position() {
        Some(position) => position.map(|axis| format!("{:.3}", f64::from(axis))),
        None => ["?".to_owned(), "?".to_owned(), "?".to_owned()],
    };

    let _ = writeln!(
        out,
        "{group} {serial} {element} {} {} {comp} {chain} {} {} {x} {y} {z} {} {} {} {auth_chain} 1",
        quote(name),
        alt(structure, atom),
        seq(residue.label_seq_id()),
        ins(residue),
        occupancy(atom),
        b_factor(atom),
        seq(residue.auth_seq_id()),
    );
}

/// A normalised chain label, unique across the structure.
///
/// The depositor's label need not be unique — a file may keep a protein and its
/// waters under the same letter and separate them another way — so writing it as
/// the normalised label would merge two chains that the file kept apart.
fn normalised_label(position: usize) -> String {
    const LETTERS: &[u8; 26] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    let mut label = String::new();
    let mut remaining = position;
    loop {
        let digit = remaining % LETTERS.len();
        let Some(letter) = LETTERS.get(digit) else {
            break;
        };
        label.insert(0, char::from(*letter));
        remaining /= LETTERS.len();
        if remaining == 0 {
            break;
        }
        remaining -= 1;
    }
    label
}

fn label(structure: &Structure, symbol: Option<pdbiox_core::symbol::SymbolId>) -> String {
    if let Some(text) = symbol.and_then(|symbol| structure.resolve(symbol))
        && !text.is_empty()
    {
        return quote(text);
    }
    ".".to_owned()
}

fn alt(structure: &Structure, atom: &AtomRef<'_>) -> String {
    let resolved = atom
        .alt_id()
        .and_then(pdbiox_core::symbol::AltId::symbol)
        .and_then(|symbol| structure.resolve(symbol));
    if let Some(text) = resolved
        && !text.is_empty()
    {
        return quote(text);
    }
    ".".to_owned()
}

fn ins(residue: &ResidueRef<'_>) -> String {
    if let Some(code) = residue.ins_code()
        && !code.is_empty()
    {
        return quote(code);
    }
    "?".to_owned()
}

fn seq(value: Option<i32>) -> String {
    match value {
        Some(value) => value.to_string(),
        None => ".".to_owned(),
    }
}

fn occupancy(atom: &AtomRef<'_>) -> String {
    match atom.occupancy() {
        Some(value) => format!("{:.2}", f64::from(value)),
        None => "?".to_owned(),
    }
}

fn b_factor(atom: &AtomRef<'_>) -> String {
    match atom.b_factor() {
        Some(value) => format!("{:.2}", f64::from(value)),
        None => "?".to_owned(),
    }
}

fn quote(text: &str) -> String {
    if needs_quoting(text) {
        format!("'{text}'")
    } else {
        text.to_owned()
    }
}

#[cfg(test)]
#[path = "write_tests.rs"]
mod tests;
