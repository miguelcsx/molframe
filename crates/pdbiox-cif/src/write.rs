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
use pdbiox_core::index::ModelIndex;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_core::topology::EntityKind;
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

    write_entry_metadata(&mut out, structure);
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

    write_entities(&mut out, structure);
    crate::write_references::write(&mut out, structure);

    out.push_str(ATOM_SITE_HEADER);
    for ordinal in 0..structure.model_count() {
        let model = ModelIndex::new(ordinal as u32);
        let Some((snapshot, local_model)) = structure.model_snapshot(model) else {
            continue;
        };
        let model_number = match snapshot
            .model(local_model)
            .and_then(pdbiox_core::structure::ModelRef::number)
        {
            Some(number) => number,
            None => match i32::try_from(ordinal.saturating_add(1)) {
                Ok(number) => number,
                Err(_) => i32::MAX,
            },
        };
        write_model(&mut out, &snapshot, local_model, model_number);
    }
    out.push_str("#\n");
    crate::write_bonds::write(&mut out, structure);
    out
}

fn write_entry_metadata(out: &mut String, structure: &Structure) {
    let entry = &structure.data().entry;
    if let Some(id) = &entry.id {
        let _ = writeln!(out, "_entry.id   {}\n#", quote(id));
    }
    if let Some(title) = &entry.title {
        let _ = writeln!(out, "_struct.title   {}\n#", quote(title));
    }
    if let Some(method) = &entry.method {
        let _ = writeln!(out, "_exptl.method   {}\n#", quote(method));
    }
    if let Some(resolution) = entry.resolution {
        let _ = writeln!(out, "_refine.ls_d_res_high   {resolution:.4}\n#");
    }
}

fn write_entities(out: &mut String, structure: &Structure) {
    let entities = &structure.data().topology.entities;
    if entities.is_empty() {
        return;
    }
    out.push_str("loop_\n_entity.id\n_entity.type\n_entity.pdbx_description\n");
    for entity in entities.iter() {
        let id = label(structure, entities.id(entity));
        let kind = match entities.kind(entity) {
            Some(EntityKind::Polymer) => "polymer",
            Some(EntityKind::NonPolymer) => "non-polymer",
            Some(EntityKind::Water) => "water",
            Some(EntityKind::Branched) => "branched",
            _ => "?",
        };
        let description = label(structure, entities.description(entity));
        let _ = writeln!(out, "{id} {kind} {description}");
    }
    out.push_str("#\n");

    let has_sequence = entities
        .iter()
        .any(|entity| !entities.canonical_sequence(entity).is_empty());
    if !has_sequence {
        return;
    }
    out.push_str(
        "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n",
    );
    for entity in entities.iter() {
        let id = label(structure, entities.id(entity));
        for (position, component) in entities.canonical_sequence(entity).iter().enumerate() {
            let component = label(structure, Some(*component));
            let _ = writeln!(out, "{id} {} {component}", position + 1);
        }
    }
    out.push_str("#\n");
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
_atom_site.label_entity_id\n\
_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
_atom_site.occupancy\n\
_atom_site.B_iso_or_equiv\n\
_atom_site.auth_seq_id\n\
_atom_site.auth_comp_id\n\
_atom_site.auth_asym_id\n\
_atom_site.auth_atom_id\n\
_atom_site.pdbx_PDB_model_num\n";

fn write_model(out: &mut String, structure: &Structure, model: ModelIndex, model_number: i32) {
    let data = structure.data();
    let mut serial = 1u32;
    for (position, chain) in data.chains().enumerate() {
        let chain_label = match chain.label() {
            Some(label) => quote(label),
            None => normalised_label(position),
        };
        let auth_label = label(
            structure,
            chain.auth_asym_id().or_else(|| chain.label_asym_id()),
        );
        let entity = chain
            .entity()
            .and_then(|entity| data.topology.entities.id(entity))
            .and_then(|symbol| structure.resolve(symbol))
            .map_or_else(|| ".".to_owned(), quote);
        let context = AtomSiteContext {
            chain: &chain_label,
            entity: &entity,
            auth_chain: &auth_label,
            model_number,
        };
        for residue in chain.residues() {
            for atom in residue.atoms() {
                let position = atom.position().and_then(|_| {
                    structure
                        .model_positions(model)
                        .and_then(|positions| positions.get(atom.index().as_usize()).copied())
                });
                write_atom(out, structure, &atom, &residue, position, &context, serial);
                serial += 1;
            }
        }
    }
}

struct AtomSiteContext<'a> {
    chain: &'a str,
    entity: &'a str,
    auth_chain: &'a str,
    model_number: i32,
}

fn write_atom(
    out: &mut String,
    structure: &Structure,
    atom: &AtomRef<'_>,
    residue: &ResidueRef<'_>,
    position: Option<[f32; 3]>,
    context: &AtomSiteContext<'_>,
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
    let comp = match atom.component_name() {
        Some(comp) => comp,
        None => "?",
    };
    let auth_comp = match residue.auth_name().or_else(|| residue.name()) {
        Some(comp) => comp,
        None => "?",
    };
    let auth_name = match atom.auth_name().or_else(|| atom.name()) {
        Some(name) => name,
        None => "?",
    };
    // An atom whose position was never recorded is written as unrecorded rather
    // than as an origin the experiment never observed.
    let [x, y, z] = match position {
        Some(position) => position.map(|axis| format!("{:.3}", f64::from(axis))),
        None => ["?".to_owned(), "?".to_owned(), "?".to_owned()],
    };

    let _ = writeln!(
        out,
        "{group} {serial} {element} {} {} {comp} {} {} {} {} {x} {y} {z} {} {} {} {} {} {} {}",
        quote(name),
        alt(structure, atom),
        context.chain,
        context.entity,
        seq(residue.label_seq_id()),
        ins(residue),
        occupancy(atom),
        b_factor(atom),
        seq(residue.auth_seq_id()),
        quote(auth_comp),
        context.auth_chain,
        quote(auth_name),
        context.model_number,
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

pub(crate) fn quote(text: &str) -> String {
    if needs_quoting(text) {
        format!("'{text}'")
    } else {
        text.to_owned()
    }
}

#[cfg(test)]
#[path = "write_tests.rs"]
mod tests;
