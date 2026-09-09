//! Assembling a structure from a document's categories.
//!
//! The coordinate category is the only one a file must have for this to produce
//! anything. Everything else — the entry's own name, the cell, the chemical
//! species — is read when present and reported as absent when not, rather than
//! being invented.

use super::atoms::{AtomBuilder, AtomSiteRowSink};
use super::ensemble::ragged_model_numbers;
use super::ragged::{RaggedBuilder, RaggedParts};
use crate::document::{CifValue, Document};
use num_traits::ToPrimitive;
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::index::EntityIndex;
use pdbiox_core::io::{ReadOptions, ReadResult};
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::structure::{CoordinateStore, Structure, StructureData, UnitCell};
use pdbiox_core::symbol::SymbolId;
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

    let ragged = ragged_model_numbers(atom_site, options.only_first_model);
    if !ragged.is_empty() {
        return Ok(lower_ragged(block, atom_site, options, &ragged));
    }

    let (structure, findings) = lower_model(block, atom_site, options, None);
    Ok((structure, findings))
}

fn lower_model(
    block: &crate::document::DataBlock,
    atom_site: &crate::document::Category,
    options: &ReadOptions,
    model: Option<i64>,
) -> (Structure, Vec<Diagnostic>) {
    let (data, findings, asym_entities) = prepare_model(block, options, Vec::new());
    let builder = AtomBuilder::new(data, findings, options, asym_entities);
    let (data, findings, coords) = match model {
        Some(model) => builder.only_model(model).read(atom_site),
        None => builder.read(atom_site),
    };
    finish_model(block, data, findings, coords)
}

pub(super) fn prepare_model(
    block: &crate::document::DataBlock,
    options: &ReadOptions,
    initial_findings: Vec<Diagnostic>,
) -> (StructureData, Diagnostics, Vec<AsymEntity>) {
    let mut findings = Diagnostics::with_capacity(initial_findings.len());
    findings.extend(initial_findings);
    let mut data = StructureData::empty();
    if !options.only_atomic_coords {
        read_entry(block, &mut data, &mut findings);
        read_cell(block, &mut data);
    }
    // Entity identity is part of topology rather than optional metadata. It is
    // needed even for a coordinate-only read so chains never become detached
    // from the chemical species declared by the file.
    read_entities(block, &mut data, &mut findings);
    let asym_entities = read_asym_entities(block, &mut data, &mut findings);
    (data, findings, asym_entities)
}

pub(super) fn finish_model(
    block: &crate::document::DataBlock,
    mut data: StructureData,
    mut findings: Diagnostics,
    coords: CoordinateStore,
) -> (Structure, Vec<Diagnostic>) {
    data.coords = coords;
    super::bonds::read(block, &mut data, &mut findings);

    for finding in pdbiox_core::structure::validate(&data) {
        findings.push(finding);
    }
    let references = super::references::read(block, &mut findings);
    let structure = Structure::new(data);
    let structure = if references.sequences.is_empty() && references.alignments.is_empty() {
        structure
    } else {
        structure.with_extension(
            pdbiox_core::structure::SEQUENCE_REFERENCES_EXTENSION,
            references,
        )
    };
    (structure, findings.finish())
}

/// Lowers borrowed coordinate rows against a metadata-only CIF document.
///
/// The caller must feed complete rows in deposition order. Every row and every
/// string it returns need remain valid only for the duration of one synchronous
/// [`AtomSiteRowSink::feed`] call. The metadata document must contain the first
/// data block and any categories needed for entry, entity, bond and reference
/// interpretation; it does not need to retain `atom_site`.
///
/// This seam is shared with binary readers so dictionary-indexed columns can be
/// lowered without allocating a [`CifValue`] for every cell.
///
/// # Errors
///
/// Returns parser or container findings from `feed`, a missing data block, or
/// findings rejected by the selected read mode.
#[doc(hidden)]
pub fn lower_atom_site_with(
    metadata: &Document,
    options: &ReadOptions,
    initial_findings: Vec<Diagnostic>,
    feed: impl FnOnce(&mut dyn AtomSiteRowSink) -> Result<(), Vec<Diagnostic>>,
) -> ReadResult {
    let Some(block) = metadata.first_block() else {
        return Err(vec![Diagnostic::new(Code::E1106)]);
    };
    let (data, findings, asym_entities) = prepare_model(block, options, initial_findings);
    let mut builder = AtomBuilder::new(data, findings, options, asym_entities);
    if let Err(errors) = feed(&mut builder) {
        let mut findings = builder.abort();
        findings.extend(errors);
        return Err(findings.finish());
    }
    let (data, findings, coords) = builder.finish();
    let (structure, findings) = finish_model(block, data, findings, coords);
    options.finish(structure, findings)
}

/// Lowers one externally classified model without retaining atom signatures.
///
/// The caller must guarantee that `feed` presents at most one selected model.
/// This contract lets the lowerer avoid an atom-sized identity side buffer.
///
/// # Errors
///
/// Returns parser or container findings from `feed`, a missing data block, or
/// findings rejected by the selected read mode.
#[doc(hidden)]
pub fn lower_single_atom_site_with(
    metadata: &Document,
    options: &ReadOptions,
    initial_findings: Vec<Diagnostic>,
    atom_capacity: usize,
    feed: impl FnOnce(&mut dyn AtomSiteRowSink) -> Result<(), Vec<Diagnostic>>,
) -> ReadResult {
    let Some(block) = metadata.first_block() else {
        return Err(vec![Diagnostic::new(Code::E1106)]);
    };
    let (data, findings, asym_entities) = prepare_model(block, options, initial_findings);
    let mut builder =
        AtomBuilder::new(data, findings, options, asym_entities).without_identity_tracking();
    builder.reserve_atoms(atom_capacity);
    if let Err(errors) = feed(&mut builder) {
        let mut findings = builder.abort();
        findings.extend(errors);
        return Err(findings.finish());
    }
    let (data, findings, coords) = builder.finish();
    let (structure, findings) = finish_model(block, data, findings, coords);
    options.finish(structure, findings)
}

/// Lowers borrowed rows whose models carry independent atom topology.
///
/// Rows must be grouped by deposited model number and fed in deposition order.
/// The sink finalises one model as the next starts, so auxiliary memory does not
/// grow with the number of source rows beyond the structures retained in the
/// ragged result itself.
///
/// # Errors
///
/// Returns parser or container findings from `feed`, a missing data block, or
/// findings rejected by the selected read mode.
#[doc(hidden)]
pub fn lower_ragged_atom_site_with(
    metadata: &Document,
    options: &ReadOptions,
    initial_findings: Vec<Diagnostic>,
    model_capacity: usize,
    feed: impl FnOnce(&mut dyn AtomSiteRowSink) -> Result<(), Vec<Diagnostic>>,
) -> ReadResult {
    let Some(block) = metadata.first_block() else {
        return Err(vec![Diagnostic::new(Code::E1106)]);
    };
    let mut findings = Diagnostics::with_capacity(initial_findings.len());
    findings.extend(initial_findings);
    let mut builder = RaggedBuilder::new(block, options, findings, model_capacity);
    if let Err(errors) = feed(&mut builder) {
        let mut findings = builder.abort();
        findings.extend(errors);
        return Err(findings.finish());
    }
    let parts = builder.finish();
    let (structure, findings) = assemble_ragged(block, options, parts);
    options.finish(structure, findings)
}

fn lower_ragged(
    block: &crate::document::DataBlock,
    atom_site: &crate::document::Category,
    options: &ReadOptions,
    model_numbers: &[i64],
) -> (Structure, Vec<Diagnostic>) {
    let mut findings = Diagnostics::new();
    let mut models = Vec::with_capacity(model_numbers.len());
    for number in model_numbers {
        let (model, model_findings) = lower_model(block, atom_site, options, Some(*number));
        findings.extend(model_findings);
        models.push(model);
    }

    assemble_ragged(
        block,
        options,
        RaggedParts {
            models,
            model_numbers: model_numbers.to_vec(),
            findings,
        },
    )
}

pub(super) fn assemble_ragged(
    block: &crate::document::DataBlock,
    options: &ReadOptions,
    parts: RaggedParts,
) -> (Structure, Vec<Diagnostic>) {
    let RaggedParts {
        models,
        model_numbers,
        mut findings,
    } = parts;
    let mut data = StructureData::empty();
    if !options.only_atomic_coords {
        read_entry(block, &mut data, &mut findings);
        read_cell(block, &mut data);
    }
    for number in model_numbers {
        let deposited = if let Ok(number) = i32::try_from(number) {
            number
        } else {
            findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("model number is outside the supported integer range")
                    .in_category(ATOM_SITE),
            );
            i32::MAX
        };
        if data.topology.models.push(deposited, 0..0).is_err() {
            findings.push(Diagnostic::new(Code::E3001));
        }
    }
    data.coords = CoordinateStore::Ragged { models };
    (Structure::new(data), findings.finish())
}

/// The entity assigned to one normalised chain identifier.
#[derive(Clone, Copy)]
pub(super) struct AsymEntity {
    pub(super) asym_id: SymbolId,
    pub(super) entity: EntityIndex,
}

/// Reads the explicit chain-to-entity relation.
pub(super) fn read_asym_entities(
    block: &crate::document::DataBlock,
    data: &mut StructureData,
    findings: &mut Diagnostics,
) -> Vec<AsymEntity> {
    let Some(asym) = block.category("struct_asym") else {
        return Vec::new();
    };
    let mut mappings = Vec::with_capacity(asym.row_count());
    for row in 0..asym.row_count() {
        let (Some(asym_id), Some(entity_id)) = (
            asym.identifier("id", row),
            asym.identifier("entity_id", row),
        ) else {
            continue;
        };
        let (Ok(asym_id), Ok(entity_id)) = (
            data.dictionary.intern(&asym_id),
            data.dictionary.intern(&entity_id),
        ) else {
            continue;
        };
        let entity = if let Some(entity) = data.topology.entities.find_by_id(entity_id) {
            entity
        } else {
            let Ok(entity) = data.topology.entities.push(
                entity_id,
                EntityKind::Unknown,
                OptionalSymbol::NONE,
                &[],
            ) else {
                findings.push(Diagnostic::new(Code::E3001));
                continue;
            };
            entity
        };
        mappings.push(AsymEntity { asym_id, entity });
    }
    mappings
}

/// Reads what the entry says about itself.
pub(super) fn read_entry(
    block: &crate::document::DataBlock,
    data: &mut StructureData,
    findings: &mut Diagnostics,
) {
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
        match resolution.to_f32() {
            Some(resolution) => data.entry.resolution = Some(resolution),
            None => findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("resolution is outside the supported floating-point range")
                    .in_category("refine"),
            ),
        }
    }
    if let Some(method) = block.category("exptl")
        && let Some(text) = method.text("method", 0)
    {
        data.entry.method = Some(text.into());
    }
}

/// Reads the crystallographic cell.
pub(super) fn read_cell(block: &crate::document::DataBlock, data: &mut StructureData) {
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
pub(super) fn read_entities(
    block: &crate::document::DataBlock,
    data: &mut StructureData,
    findings: &mut Diagnostics,
) {
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
        if data
            .topology
            .entities
            .push(id, kind, description, &sequence)
            .is_err()
        {
            findings.push(Diagnostic::new(Code::E3001));
            return;
        }
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
