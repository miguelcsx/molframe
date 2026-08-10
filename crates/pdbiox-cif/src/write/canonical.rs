//! Fallible canonical projection from semantic structure data to mmCIF.
//!
//! Validation completes before output allocation, so failure never returns a
//! partial file. Values that CIF explicitly permits to be unknown remain
//! sentinels; identifiers required to address atoms, models, and bond endpoints
//! are never invented implicitly.

use super::options::{CifWriteError, CifWriteOptions, valid_block_id};
use super::value::quote_text;
use pdbiox_core::index::ModelIndex;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_core::topology::EntityKind;
use std::fmt::Write as _;

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
    let block_id = preflight(structure, options)?;
    let mut out = String::with_capacity(structure.atom_count() as usize * 100);
    let _ = writeln!(out, "data_{block_id}");
    out.push_str("#\n");
    write_entry_metadata(&mut out, structure);
    write_cell(&mut out, structure);
    write_entities(&mut out, structure);
    super::references::write(&mut out, structure);
    write_atoms(&mut out, structure)?;
    super::bonds::write(&mut out, structure, options)?;
    Ok(out)
}

fn preflight<'a>(
    structure: &'a Structure,
    options: &'a CifWriteOptions,
) -> Result<&'a str, CifWriteError> {
    let Some(block_id) = options.block_id().or(structure.data().entry.id.as_deref()) else {
        return Err(CifWriteError::MissingBlockId);
    };
    if !valid_block_id(block_id) {
        return Err(CifWriteError::InvalidBlockId(block_id.to_owned()));
    }
    for position in 0..structure.model_count() {
        let index = model_index(position)?;
        if structure
            .model(index)
            .and_then(pdbiox_core::structure::ModelRef::number)
            .is_none()
        {
            return Err(CifWriteError::MissingModelNumber { model: index.get() });
        }
    }
    for position in 0..structure.atom_count() {
        let Some(atom) = structure.atom(pdbiox_core::index::AtomIndex::new(position)) else {
            return Err(CifWriteError::MissingAtomField {
                atom: position,
                field: "atom row",
            });
        };
        validate_atom(structure, atom)?;
    }
    super::bonds::preflight(structure, options)?;
    Ok(block_id)
}

fn validate_atom(structure: &Structure, atom: AtomRef<'_>) -> Result<(), CifWriteError> {
    let atom_index = atom.index().get();
    if atom.atom_site_id().is_none_or(|id| id == 0) {
        return Err(CifWriteError::MissingAtomSiteId { atom: atom_index });
    }
    required(atom.name(), atom_index, "label_atom_id")?;
    required(atom.component_name(), atom_index, "label_comp_id")?;
    let chain = atom
        .residue()
        .and_then(|residue| chain_of(structure, residue))
        .and_then(pdbiox_core::structure::ChainRef::label);
    required(chain, atom_index, "label_asym_id")?;
    Ok(())
}

fn required(value: Option<&str>, atom: u32, field: &'static str) -> Result<(), CifWriteError> {
    if value.is_some_and(|text| !text.is_empty()) {
        Ok(())
    } else {
        Err(CifWriteError::MissingAtomField { atom, field })
    }
}

fn write_entry_metadata(out: &mut String, structure: &Structure) {
    let entry = &structure.data().entry;
    if let Some(id) = &entry.id {
        let _ = writeln!(out, "_entry.id   {}\n#", quote_text(id));
    }
    if let Some(title) = &entry.title {
        let _ = writeln!(out, "_struct.title   {}\n#", quote_text(title));
    }
    if let Some(method) = &entry.method {
        let _ = writeln!(out, "_exptl.method   {}\n#", quote_text(method));
    }
    if let Some(resolution) = entry.resolution {
        let _ = writeln!(out, "_refine.ls_d_res_high   {resolution:.4}\n#");
    }
}

fn write_cell(out: &mut String, structure: &Structure) {
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

fn write_entities(out: &mut String, structure: &Structure) {
    let entities = &structure.data().topology.entities;
    if entities.is_empty() {
        return;
    }
    out.push_str("loop_\n_entity.id\n_entity.type\n_entity.pdbx_description\n");
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
    out.push_str("#\n");
    write_entity_sequences(out, structure);
}

fn write_entity_sequences(out: &mut String, structure: &Structure) {
    let entities = &structure.data().topology.entities;
    if !entities
        .iter()
        .any(|entity| !entities.canonical_sequence(entity).is_empty())
    {
        return;
    }
    out.push_str(
        "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n",
    );
    for entity in entities.iter() {
        let id = symbol_or_dot(structure, entities.id(entity));
        for (position, component) in entities.canonical_sequence(entity).iter().enumerate() {
            let component = symbol_or_dot(structure, Some(*component));
            let _ = writeln!(out, "{id} {} {component}", position + 1);
        }
    }
    out.push_str("#\n");
}

const ATOM_SITE_HEADER: &str = "loop_\n\
_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n\
_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n\
_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n";

fn write_atoms(out: &mut String, structure: &Structure) -> Result<(), CifWriteError> {
    out.push_str(ATOM_SITE_HEADER);
    for position in 0..structure.model_count() {
        let model = model_index(position)?;
        let Some((snapshot, local_model)) = structure.model_snapshot(model) else {
            return Err(CifWriteError::MissingModelNumber { model: model.get() });
        };
        let Some(model_number) = snapshot
            .model(local_model)
            .and_then(pdbiox_core::structure::ModelRef::number)
        else {
            return Err(CifWriteError::MissingModelNumber { model: model.get() });
        };
        write_model(out, &snapshot, local_model, model_number)?;
    }
    out.push_str("#\n");
    Ok(())
}

fn model_index(position: usize) -> Result<ModelIndex, CifWriteError> {
    u32::try_from(position)
        .map(ModelIndex::new)
        .map_err(|_| CifWriteError::ModelIndexOverflow { model: position })
}

fn write_model(
    out: &mut String,
    structure: &Structure,
    model: ModelIndex,
    model_number: i32,
) -> Result<(), CifWriteError> {
    for chain in structure.data().chains() {
        let chain_label = chain.label().map(quote_text);
        let auth_label = optional_unknown(chain.auth_label());
        let entity = chain
            .entity()
            .and_then(|entity| structure.data().topology.entities.id(entity));
        let entity = symbol_or_dot(structure, entity);
        let context = AtomSiteContext {
            chain: chain_label.as_deref(),
            entity: &entity,
            auth_chain: &auth_label,
            model_number,
        };
        for residue in chain.residues() {
            for atom in residue.atoms() {
                let position = structure
                    .model_positions(model)
                    .and_then(|positions| positions.get(atom.index().as_usize()).copied());
                write_atom(out, structure, atom, residue, position, &context)?;
            }
        }
    }
    Ok(())
}

struct AtomSiteContext<'a> {
    chain: Option<&'a str>,
    entity: &'a str,
    auth_chain: &'a str,
    model_number: i32,
}

fn write_atom(
    out: &mut String,
    structure: &Structure,
    atom: AtomRef<'_>,
    residue: ResidueRef<'_>,
    position: Option<[f32; 3]>,
    context: &AtomSiteContext<'_>,
) -> Result<(), CifWriteError> {
    let atom_index = atom.index().get();
    let Some(atom_id) = atom.atom_site_id() else {
        return Err(CifWriteError::MissingAtomSiteId { atom: atom_index });
    };
    let Some(atom_name) = atom.name() else {
        return Err(CifWriteError::MissingAtomField {
            atom: atom_index,
            field: "label_atom_id",
        });
    };
    let Some(component_name) = atom.component_name() else {
        return Err(CifWriteError::MissingAtomField {
            atom: atom_index,
            field: "label_comp_id",
        });
    };
    let Some(chain_label) = context.chain else {
        return Err(CifWriteError::MissingAtomField {
            atom: atom_index,
            field: "label_asym_id",
        });
    };
    let group = if residue.is_het() { "HETATM" } else { "ATOM" };
    let element = atom
        .element()
        .map_or_else(|| "?".to_owned(), |item| item.symbol().to_uppercase());
    let [x, y, z] = position.map_or_else(
        || ["?".to_owned(), "?".to_owned(), "?".to_owned()],
        |point| point.map(|axis| format!("{:.3}", f64::from(axis))),
    );
    let _ = writeln!(
        out,
        "{group} {} {element} {} {} {} {} {} {} {} {x} {y} {z} {} {} {} {} {} {} {}",
        atom_id,
        quote_text(atom_name),
        alternate(structure, atom),
        quote_text(component_name),
        chain_label,
        context.entity,
        sequence(residue.label_seq_id()),
        insertion(residue),
        optional_number(atom.occupancy()),
        optional_number(atom.b_factor()),
        sequence(residue.auth_seq_id()),
        optional_unknown(residue.auth_name()),
        context.auth_chain,
        optional_unknown(atom.auth_name()),
        context.model_number,
    );
    Ok(())
}

fn symbol_or_dot(structure: &Structure, symbol: Option<pdbiox_core::symbol::SymbolId>) -> String {
    match symbol.and_then(|id| structure.resolve(id)) {
        Some(text) if !text.is_empty() => quote_text(text),
        Some(_) | None => ".".to_owned(),
    }
}

fn optional_unknown(value: Option<&str>) -> String {
    match value {
        Some(text) if !text.is_empty() => quote_text(text),
        Some(_) | None => "?".to_owned(),
    }
}

fn alternate(structure: &Structure, atom: AtomRef<'_>) -> String {
    let value = atom
        .alt_id()
        .and_then(pdbiox_core::symbol::AltId::symbol)
        .and_then(|id| structure.resolve(id));
    match value {
        Some(text) if !text.is_empty() => quote_text(text),
        Some(_) | None => ".".to_owned(),
    }
}

fn insertion(residue: ResidueRef<'_>) -> String {
    match residue.ins_code() {
        Some(code) if !code.is_empty() => quote_text(code),
        Some(_) | None => "?".to_owned(),
    }
}

fn sequence(value: Option<i32>) -> String {
    match value {
        Some(value) => value.to_string(),
        None => ".".to_owned(),
    }
}

fn optional_number(value: Option<f32>) -> String {
    match value {
        Some(value) => format!("{value:.2}"),
        None => "?".to_owned(),
    }
}

fn chain_of<'a>(
    structure: &'a Structure,
    residue: ResidueRef<'_>,
) -> Option<pdbiox_core::structure::ChainRef<'a>> {
    let index = structure
        .data()
        .topology
        .chains
        .containing(residue.index().get())?;
    structure.chain(index)
}

#[cfg(test)]
#[path = "canonical_tests.rs"]
mod tests;
