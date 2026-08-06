//! Assembling a structure from a document's categories.
//!
//! The coordinate category is the only one a file must have for this to produce
//! anything. Everything else — the entry's own name, the cell, the chemical
//! species — is read when present and reported as absent when not, rather than
//! being invented.

use super::atoms::AtomBuilder;
use crate::document::{CifValue, Document};
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::io::{ReadOptions, ReadResult};
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::structure::{Structure, StructureData, UnitCell};
use pdbiox_core::topology::EntityKind;

/// The category coordinates live in.
const ATOM_SITE: &str = "atom_site";

/// Builds a structure from the first block of a document.
///
/// # Errors
///
/// Returns the findings that stopped it: a document with no block, or a block
/// with no coordinates.
pub fn lower(document: &Document, options: &ReadOptions) -> ReadResult {
    let mut findings = Diagnostics::new();
    let Some(block) = document.first_block() else {
        return Err(vec![Diagnostic::new(Code::E1106)]);
    };
    let Some(atom_site) = block.category(ATOM_SITE) else {
        return Err(vec![
            Diagnostic::new(Code::E2001)
                .with_message("the block contains no coordinates")
                .in_category(ATOM_SITE),
        ]);
    };

    let mut data = StructureData::empty();
    if !options.only_atomic_coords {
        read_entry(block, &mut data);
        read_cell(block, &mut data);
        read_entities(block, &mut data);
    }

    let coords = AtomBuilder::new(&mut data, &mut findings, options).read(atom_site);
    data.coords = coords;

    for finding in pdbiox_core::structure::validate(&data) {
        findings.push(finding);
    }
    Ok((Structure::new(data), findings.finish()))
}

/// Reads what the entry says about itself.
fn read_entry(block: &crate::document::DataBlock, data: &mut StructureData) {
    if let Some(entry) = block.category("entry")
        && let Some(id) = entry.text("id", 0)
    {
        data.entry.id = Some(id.into());
    } else {
        // A block name is the identifier a file without an entry category has.
        let name = block.name();
        if !name.is_empty() {
            data.entry.id = Some(name.into());
        }
    }
    if let Some(struct_category) = block.category("struct")
        && let Some(title) = struct_category.text("title", 0)
    {
        data.entry.title = Some(title.into());
    }
    if let Some(refine) = block.category("refine")
        && let Some(resolution) = refine
            .value("ls_d_res_high", 0)
            .and_then(CifValue::as_float)
    {
        data.entry.resolution = Some(resolution as f32);
    }
    if let Some(method) = block.category("exptl")
        && let Some(text) = method.text("method", 0)
    {
        data.entry.method = Some(text.into());
    }
}

/// Reads the crystallographic cell.
fn read_cell(block: &crate::document::DataBlock, data: &mut StructureData) {
    let Some(cell) = block.category("cell") else {
        return;
    };
    let read = |item: &str| cell.value(item, 0).and_then(CifValue::as_float);
    let (Some(a), Some(b), Some(c)) = (read("length_a"), read("length_b"), read("length_c")) else {
        return;
    };
    let (Some(alpha), Some(beta), Some(gamma)) =
        (read("angle_alpha"), read("angle_beta"), read("angle_gamma"))
    else {
        return;
    };
    data.cell = Some(UnitCell {
        lengths: [a, b, c],
        angles: [alpha, beta, gamma],
    });
}

/// Reads the distinct chemical species, and the sequence each should have.
fn read_entities(block: &crate::document::DataBlock, data: &mut StructureData) {
    let Some(entities) = block.category("entity") else {
        return;
    };
    for row in 0..entities.row_count() {
        let Some(id) = entities.identifier("id", row) else {
            continue;
        };
        let kind = match entities.text("type", row) {
            Some("polymer") => EntityKind::Polymer,
            Some("non-polymer") => EntityKind::NonPolymer,
            Some("water") => EntityKind::Water,
            Some("branched") => EntityKind::Branched,
            _ => EntityKind::Unknown,
        };
        let Ok(id) = data.dictionary.intern(&id) else {
            continue;
        };
        let description = match entities.text("pdbx_description", row) {
            Some(text) => match data.dictionary.intern(text) {
                Ok(symbol) => OptionalSymbol::some(symbol),
                Err(_) => OptionalSymbol::NONE,
            },
            None => OptionalSymbol::NONE,
        };
        let sequence = canonical_sequence(block, data, entities.identifier("id", row).as_deref());
        data.topology
            .entities
            .push(id, kind, description, &sequence);
    }
}

/// The sequence an entity should have, from the category that lists it.
fn canonical_sequence(
    block: &crate::document::DataBlock,
    data: &mut StructureData,
    entity_id: Option<&str>,
) -> Vec<pdbiox_core::symbol::SymbolId> {
    let Some(entity_id) = entity_id else {
        return Vec::new();
    };
    let Some(sequence) = block.category("entity_poly_seq") else {
        return Vec::new();
    };

    let mut components = Vec::new();
    for row in 0..sequence.row_count() {
        if sequence.identifier("entity_id", row).as_deref() != Some(entity_id) {
            continue;
        }
        let Some(name) = sequence.identifier("mon_id", row) else {
            continue;
        };
        if let Ok(symbol) = data.dictionary.intern(&name) {
            components.push(symbol);
        }
    }
    components
}
