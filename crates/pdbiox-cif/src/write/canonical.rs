//! Fallible canonical projection from semantic structure data to mmCIF.
//!
//! Validation completes before output allocation, so failure never returns a
//! partial file. Values that CIF explicitly permits to be unknown remain
//! sentinels; identifiers required to address atoms, models, and bond endpoints
//! are never invented implicitly.

use super::options::{CifWriteError, CifWriteOptions, CifWriteToError};
use super::projection::{
    CanonicalAtomRow, CanonicalProjection, CanonicalValue, canonical_projection,
};
use super::value::quoted;
use pdbiox_core::io::TextOutput;
use pdbiox_core::structure::Structure;
use pdbiox_core::topology::EntityKind;
use std::fmt::{self, Display, Formatter, Write as _};
use std::io::Write as IoWrite;

/// Writes a structure using only identifiers retained in the structure.
///
/// # Errors
///
/// Returns [`CifWriteError`] before allocating output when a required identity
/// is absent or connectivity would require generated connection identifiers.
pub fn write_canonical(structure: &Structure) -> Result<String, CifWriteError> {
    write_canonical_with_options(structure, &CifWriteOptions::new())
}

/// Writes a structure with explicit choices for identities not retained by the graph.
///
/// # Errors
///
/// Returns [`CifWriteError`] before allocating output when preflight finds an
/// unrepresentable atom, model, bond endpoint, or block identifier.
pub fn write_canonical_with_options(
    structure: &Structure,
    options: &CifWriteOptions,
) -> Result<String, CifWriteError> {
    let mut out = String::with_capacity(structure.atom_count() as usize * 100);
    render_canonical(&mut out, structure, options)?;
    Ok(out)
}

/// Streams canonical mmCIF directly to a byte destination.
///
/// Projection validation completes before the first write. The writer retains
/// no atom rows and uses only formatting-sized temporary storage.
///
/// # Errors
///
/// Returns a canonical projection refusal or the destination I/O error.
pub fn write_canonical_to<W: IoWrite>(
    structure: &Structure,
    options: &CifWriteOptions,
    output: &mut W,
) -> Result<(), CifWriteToError> {
    let mut output = TextOutput::new(output);
    render_canonical(&mut output, structure, options)?;
    output.finish().map_err(CifWriteToError::Output)
}

fn render_canonical(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &CifWriteOptions,
) -> Result<(), CifWriteError> {
    let projection = canonical_projection(structure, options)?;
    let _ = writeln!(out, "data_{}", projection.block_id());
    let _ = out.write_str("#\n");
    write_entry_metadata(out, structure);
    write_cell(out, structure);
    write_entities(out, structure);
    super::references::write(out, structure);
    write_atoms(out, projection)?;
    super::bonds::write(out, structure, options)?;
    Ok(())
}

fn write_entry_metadata(out: &mut impl fmt::Write, structure: &Structure) {
    let entry = &structure.data().entry;
    if let Some(id) = &entry.id {
        let _ = writeln!(out, "_entry.id   {}\n#", quoted(id));
    }
    if let Some(title) = &entry.title {
        let _ = writeln!(out, "_struct.title   {}\n#", quoted(title));
    }
    if let Some(method) = &entry.method {
        let _ = writeln!(out, "_exptl.method   {}\n#", quoted(method));
    }
    if let Some(resolution) = entry.resolution {
        let _ = writeln!(out, "_refine.ls_d_res_high   {resolution:.4}\n#");
    }
}

fn write_cell(out: &mut impl fmt::Write, structure: &Structure) {
    let Some(cell) = structure.data().cell else {
        return;
    };
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

fn write_entities(out: &mut impl fmt::Write, structure: &Structure) {
    let entities = &structure.data().topology.entities;
    if entities.is_empty() {
        return;
    }
    let _ = out.write_str("loop_\n_entity.id\n_entity.type\n_entity.pdbx_description\n");
    for entity in entities.iter() {
        let id = symbol_or_dot(structure, entities.id(entity));
        let kind = match entities.kind(entity) {
            Some(EntityKind::Polymer) => "polymer",
            Some(EntityKind::NonPolymer) => "non-polymer",
            Some(EntityKind::Water) => "water",
            Some(EntityKind::Branched) => "branched",
            _ => "?",
        };
        let description = symbol_or_dot(structure, entities.description(entity));
        let _ = writeln!(out, "{id} {kind} {description}");
    }
    let _ = out.write_str("#\n");
    write_entity_sequences(out, structure);
}

fn write_entity_sequences(out: &mut impl fmt::Write, structure: &Structure) {
    let entities = &structure.data().topology.entities;
    if !entities
        .iter()
        .any(|entity| !entities.canonical_sequence(entity).is_empty())
    {
        return;
    }
    let _ = out.write_str(
        "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n",
    );
    for entity in entities.iter() {
        let id = symbol_or_dot(structure, entities.id(entity));
        for (position, component) in entities.canonical_sequence(entity).iter().enumerate() {
            let component = symbol_or_dot(structure, Some(*component));
            let _ = writeln!(out, "{id} {} {component}", position + 1);
        }
    }
    let _ = out.write_str("#\n");
}

const ATOM_SITE_HEADER: &str = "loop_\n\
_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n\
_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n\
_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n";

fn write_atoms(
    out: &mut impl fmt::Write,
    projection: CanonicalProjection<'_>,
) -> Result<(), CifWriteError> {
    let _ = out.write_str(ATOM_SITE_HEADER);
    projection.visit_atom_rows(|row| {
        let _ = write_atom(out, row);
    })?;
    let _ = out.write_str("#\n");
    Ok(())
}

fn write_atom(out: &mut impl fmt::Write, row: CanonicalAtomRow<'_>) -> fmt::Result {
    let element = match row.element() {
        CanonicalValue::Present(item) => ElementValue::Present(item.symbol()),
        CanonicalValue::Inapplicable => ElementValue::Inapplicable,
        CanonicalValue::Unknown => ElementValue::Unknown,
    };
    writeln!(
        out,
        "{} {} {element} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        row.group(),
        row.atom_site_id(),
        quoted(row.label_atom_id()),
        text_value(row.label_alt_id()),
        quoted(row.label_comp_id()),
        quoted(row.label_asym_id()),
        text_value(row.label_entity_id()),
        integer_value(row.label_seq_id()),
        text_value(row.insertion_code()),
        float_value(row.coordinate(0), 3),
        float_value(row.coordinate(1), 3),
        float_value(row.coordinate(2), 3),
        float_value(row.occupancy(), 2),
        float_value(row.b_factor(), 2),
        integer_value(row.auth_seq_id()),
        text_value(row.auth_comp_id()),
        text_value(row.auth_asym_id()),
        text_value(row.auth_atom_id()),
        row.model_number(),
    )
}

fn text_value(value: CanonicalValue<&str>) -> TextValue<'_> {
    TextValue(value)
}

fn integer_value(value: CanonicalValue<i64>) -> IntegerValue {
    IntegerValue(value)
}

fn float_value(value: CanonicalValue<f64>, precision: usize) -> FloatValue {
    FloatValue { value, precision }
}

fn symbol_or_dot(
    structure: &Structure,
    symbol: Option<pdbiox_core::symbol::SymbolId>,
) -> TextValue<'_> {
    match symbol.and_then(|id| structure.resolve(id)) {
        Some(text) if !text.is_empty() => TextValue(CanonicalValue::Present(text)),
        Some(_) | None => TextValue(CanonicalValue::Inapplicable),
    }
}

struct TextValue<'a>(CanonicalValue<&'a str>);

impl Display for TextValue<'_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self.0 {
            CanonicalValue::Present(text) => Display::fmt(&quoted(text), formatter),
            CanonicalValue::Inapplicable => formatter.write_str("."),
            CanonicalValue::Unknown => formatter.write_str("?"),
        }
    }
}

struct IntegerValue(CanonicalValue<i64>);

impl Display for IntegerValue {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self.0 {
            CanonicalValue::Present(number) => Display::fmt(&number, formatter),
            CanonicalValue::Inapplicable => formatter.write_str("."),
            CanonicalValue::Unknown => formatter.write_str("?"),
        }
    }
}

struct FloatValue {
    value: CanonicalValue<f64>,
    precision: usize,
}

impl Display for FloatValue {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self.value {
            CanonicalValue::Present(number) => {
                write!(
                    formatter,
                    "{number:.precision$}",
                    precision = self.precision
                )
            }
            CanonicalValue::Inapplicable => formatter.write_str("."),
            CanonicalValue::Unknown => formatter.write_str("?"),
        }
    }
}

enum ElementValue<'a> {
    Present(&'a str),
    Inapplicable,
    Unknown,
}

impl Display for ElementValue<'_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Present(symbol) => {
                for byte in symbol.bytes() {
                    formatter.write_char(char::from(byte.to_ascii_uppercase()))?;
                }
                Ok(())
            }
            Self::Inapplicable => formatter.write_str("."),
            Self::Unknown => formatter.write_str("?"),
        }
    }
}

#[cfg(test)]
#[path = "canonical_tests.rs"]
mod tests;
