use super::codec::{decode_chars, decode_f32, decode_i32, decode_strings};
use super::schema::{Entity, File, Group};
use super::{
    MMTF_METADATA_EXTENSION, MmtfEntityMetadata, MmtfGroupMetadata, MmtfMetadata, MmtfOptionalField,
};
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::column::Presence;
use pdbiox_core::coords::CoordinateBlock;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::index::{AtomIndex, EntityIndex, ResidueIndex};
use pdbiox_core::io::{InputBuffer, ReadOptions, ReadResult};
use pdbiox_core::optional::{OptionalI32, OptionalSymbol};
use pdbiox_core::structure::{CoordinateStore, Structure, StructureData, UnitCell};
use pdbiox_core::symbol::{AltId, SymbolId};
use pdbiox_core::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};
use pdbiox_core::{BondOrder, BondProvenance, BondRecord, BondTableBuilder, Element};

/// Reads an MMTF `MessagePack` structure through all standard codecs.
///
/// # Errors
///
/// Returns diagnostics for malformed `MessagePack`, unsupported major versions,
/// invalid codecs, inconsistent declared counts, or model topology changes.
pub fn read_mmtf(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    let file: File = match rmp_serde::from_slice(input.as_bytes()) {
        Ok(file) => file,
        Err(error) => {
            return Err(vec![
                Diagnostic::new(Code::E1102)
                    .with_message("MMTF MessagePack could not be decoded")
                    .with_context("decoder", error.to_string()),
            ]);
        }
    };
    match lower(file, options) {
        Ok(structure) => options.finish(structure, []),
        Err(finding) => Err(vec![finding]),
    }
}

struct Decoded {
    file: File,
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    b_factor: Option<Vec<f32>>,
    occupancy: Option<Vec<f32>>,
    atom_id: Option<Vec<i32>>,
    alt_loc: Option<Vec<u8>>,
    group_id: Vec<i32>,
    group_type: Vec<i32>,
    ins_code: Option<Vec<u8>>,
    sequence_index: Option<Vec<i32>>,
    chain_id: Vec<String>,
    chain_name: Option<Vec<String>>,
    global_bonds: Vec<i32>,
    global_orders: Vec<i32>,
    global_resonance: Vec<i32>,
}

fn lower(file: File, options: &ReadOptions) -> Result<Structure, Diagnostic> {
    validate_version(&file.mmtf_version)?;
    validate_declared_counts(&file, options)?;
    let metadata = metadata(&file);
    let decoded = decode(file)?;
    validate_decoded_lengths(&decoded)?;
    build(&decoded, options)
        .map(|structure| structure.with_extension(MMTF_METADATA_EXTENSION, metadata))
}

fn metadata(file: &File) -> MmtfMetadata {
    MmtfMetadata {
        space_group: file.space_group.as_deref().map(Into::into),
        groups: file
            .group_list
            .iter()
            .map(|group| MmtfGroupMetadata {
                name: group.name.as_str().into(),
                atom_names: group
                    .atom_name_list
                    .iter()
                    .map(|name| name.as_str().into())
                    .collect(),
                elements: group.element_list.as_ref().map(|elements| {
                    elements
                        .iter()
                        .map(|element| element.as_str().into())
                        .collect()
                }),
                single_letter_code: group.single_letter_code.as_str().into(),
                chem_comp_type: group.chem_comp_type.as_str().into(),
            })
            .collect(),
        entities: file
            .entity_list
            .as_deref()
            .map_or(&[][..], |entities| entities)
            .iter()
            .map(|entity| MmtfEntityMetadata {
                description: entity.description.as_str().into(),
                kind: entity.kind.as_str().into(),
                sequence: entity.sequence.as_str().into(),
            })
            .collect(),
        optional_fields: optional_fields(file),
    }
}

fn optional_fields(file: &File) -> std::collections::BTreeSet<MmtfOptionalField> {
    [
        (file.b_factor_list.is_some(), MmtfOptionalField::BFactor),
        (file.occupancy_list.is_some(), MmtfOptionalField::Occupancy),
        (file.atom_id_list.is_some(), MmtfOptionalField::AtomId),
        (file.alt_loc_list.is_some(), MmtfOptionalField::AltLoc),
        (file.ins_code_list.is_some(), MmtfOptionalField::InsCode),
        (
            file.sequence_index_list.is_some(),
            MmtfOptionalField::SequenceIndex,
        ),
        (file.chain_name_list.is_some(), MmtfOptionalField::ChainName),
        (file.entity_list.is_some(), MmtfOptionalField::EntityList),
    ]
    .into_iter()
    .filter_map(|(present, field)| present.then_some(field))
    .collect()
}

fn decode(file: File) -> Result<Decoded, Diagnostic> {
    let x = decode_f32(&file.x_coord_list)?;
    let y = decode_f32(&file.y_coord_list)?;
    let z = decode_f32(&file.z_coord_list)?;
    let b_factor = file
        .b_factor_list
        .as_ref()
        .map(|bytes| decode_f32(bytes.as_ref()))
        .transpose()?;
    let occupancy = file
        .occupancy_list
        .as_ref()
        .map(|bytes| decode_f32(bytes.as_ref()))
        .transpose()?;
    let atom_id = file
        .atom_id_list
        .as_ref()
        .map(|bytes| decode_i32(bytes.as_ref()))
        .transpose()?;
    let alt_loc = file
        .alt_loc_list
        .as_ref()
        .map(|bytes| decode_chars(bytes.as_ref()))
        .transpose()?;
    let group_id = decode_i32(&file.group_id_list)?;
    let group_type = decode_i32(&file.group_type_list)?;
    let ins_code = file
        .ins_code_list
        .as_ref()
        .map(|bytes| decode_chars(bytes.as_ref()))
        .transpose()?;
    let sequence_index = file
        .sequence_index_list
        .as_ref()
        .map(|bytes| decode_i32(bytes.as_ref()))
        .transpose()?;
    let chain_id = decode_strings(&file.chain_id_list)?;
    let chain_name = file
        .chain_name_list
        .as_ref()
        .map(|bytes| decode_strings(bytes.as_ref()))
        .transpose()?;
    let global_bonds = decode_optional_i32(file.bond_atom_list.as_deref().map(Vec::as_slice))?;
    let global_orders = decode_optional_i32(file.bond_order_list.as_deref().map(Vec::as_slice))?;
    let global_resonance =
        decode_optional_i32(file.bond_resonance_list.as_deref().map(Vec::as_slice))?;
    Ok(Decoded {
        file,
        x,
        y,
        z,
        b_factor,
        occupancy,
        atom_id,
        alt_loc,
        group_id,
        group_type,
        ins_code,
        sequence_index,
        chain_id,
        chain_name,
        global_bonds,
        global_orders,
        global_resonance,
    })
}

fn decode_optional_i32(bytes: Option<&[u8]>) -> Result<Vec<i32>, Diagnostic> {
    match bytes {
        Some(bytes) => decode_i32(bytes),
        None => Ok(Vec::new()),
    }
}

fn build(decoded: &Decoded, options: &ReadOptions) -> Result<Structure, Diagnostic> {
    let ranges = model_ranges(decoded)?;
    let requested = if options.only_first_model {
        1
    } else {
        ranges.len()
    };
    let mut models = Vec::with_capacity(requested);
    for (model, range) in ranges.into_iter().take(requested).enumerate() {
        models.push(build_model(decoded, options, model, range)?);
    }
    if models.len() == 1 {
        return models
            .pop()
            .ok_or_else(|| schema_error("MMTF model was not built"));
    }
    if models
        .iter()
        .skip(1)
        .all(|model| same_topology(&models[0], model))
    {
        dense_ensemble(&models)
    } else {
        ragged_ensemble(models)
    }
}

#[derive(Clone, Copy)]
struct ModelRange {
    chain_start: usize,
    chain_count: usize,
    group_start: usize,
    group_count: usize,
    atom_start: usize,
    atom_count: usize,
}

fn model_ranges(decoded: &Decoded) -> Result<Vec<ModelRange>, Diagnostic> {
    let mut ranges = Vec::with_capacity(decoded.file.chains_per_model.len());
    let mut chain_start = 0usize;
    let mut group_start = 0usize;
    let mut atom_start = 0usize;
    for chain_count in &decoded.file.chains_per_model {
        let chain_count = usize_of(*chain_count)?;
        let chain_end = chain_start
            .checked_add(chain_count)
            .ok_or_else(|| schema_error("MMTF chain range overflows"))?;
        let group_count = decoded
            .file
            .groups_per_chain
            .get(chain_start..chain_end)
            .ok_or_else(|| schema_error("chainsPerModel exceeds groupsPerChain"))?
            .iter()
            .try_fold(0usize, |sum, value| {
                sum.checked_add(usize_of(*value)?)
                    .ok_or_else(|| schema_error("MMTF group count overflows"))
            })?;
        let atom_count = atoms_in_groups(decoded, group_start, group_count)?;
        ranges.push(ModelRange {
            chain_start,
            chain_count,
            group_start,
            group_count,
            atom_start,
            atom_count,
        });
        chain_start = chain_end;
        group_start += group_count;
        atom_start += atom_count;
    }
    if chain_start != usize_of(decoded.file.num_chains)?
        || group_start != usize_of(decoded.file.num_groups)?
        || atom_start != usize_of(decoded.file.num_atoms)?
    {
        return Err(schema_error(
            "MMTF hierarchy counts do not sum to declared totals",
        ));
    }
    Ok(ranges)
}

fn build_model(
    decoded: &Decoded,
    options: &ReadOptions,
    model: usize,
    range: ModelRange,
) -> Result<Structure, Diagnostic> {
    let mut data = StructureData::empty();
    data.entry.id = decoded.file.structure_id.clone().map(Into::into);
    data.entry.title = decoded.file.title.clone().map(Into::into);
    data.entry.method = decoded
        .file
        .experimental_methods
        .as_ref()
        .map(|methods| methods.join(", ").into_boxed_str());
    data.entry.resolution = decoded.file.resolution;
    data.cell = unit_cell(decoded.file.unit_cell.as_deref())?;
    let entity_by_chain = build_entities(&mut data, decoded)?;
    let mut builder = ChunkBuilder::new();
    let mut keep_map = vec![None; usize_of(decoded.file.num_atoms)?];
    let mut bond_builder = BondTableBuilder::new();
    let mut group_cursor = range.group_start;
    let mut atom_cursor = range.atom_start;
    let mut kept_atoms = 0u32;
    for chain_position in range.chain_start..range.chain_start + range.chain_count {
        let first_residue = u32_of(data.topology.residues.len(), "residue index")?;
        let group_count = usize_of(decoded.file.groups_per_chain[chain_position])?;
        for _ in 0..group_count {
            let group = group_at(decoded, group_cursor)?;
            let first_kept = kept_atoms;
            let residue_index =
                ResidueIndex::new(u32_of(data.topology.residues.len(), "residue index")?);
            for local in 0..group.atom_name_list.len() {
                let original = atom_cursor + local;
                let element = group_element(group, local);
                if options.discard_hydrogens && element.is_hydrogen() {
                    continue;
                }
                keep_map[original] = Some(AtomIndex::new(kept_atoms));
                builder.push(atom_record(
                    &mut data,
                    decoded,
                    group,
                    local,
                    original,
                    residue_index,
                )?);
                kept_atoms += 1;
            }
            let record = residue_record(&mut data, decoded, group, group_cursor, chain_position)?;
            data.topology
                .residues
                .push(record, first_kept..kept_atoms)
                .map_err(|error| {
                    schema_error("residue table rejected a row")
                        .with_context("cause", error.to_string())
                })?;
            add_group_bonds(group, atom_cursor, &keep_map, &mut bond_builder)?;
            atom_cursor += group.atom_name_list.len();
            group_cursor += 1;
        }
        let entity = *entity_by_chain
            .get(chain_position)
            .ok_or_else(|| schema_error("chain has no entity mapping"))?;
        let kind = match data.topology.entities.kind(entity) {
            Some(kind) => kind,
            None => EntityKind::Unknown,
        };
        let label = intern(&mut data, &decoded.chain_id[chain_position])?;
        let auth = decoded
            .chain_name
            .as_ref()
            .and_then(|names| names.get(chain_position))
            .map(|name| intern(&mut data, name))
            .transpose()?
            .map_or(OptionalSymbol::NONE, OptionalSymbol::some);
        data.topology
            .chains
            .push(
                ChainRecord {
                    label_asym_id: label,
                    auth_asym_id: auth,
                    entity,
                    polymer_kind: if kind == EntityKind::Polymer {
                        PolymerKind::Other
                    } else {
                        PolymerKind::None
                    },
                },
                first_residue..u32_of(data.topology.residues.len(), "residue count")?,
            )
            .map_err(|error| {
                schema_error("chain table rejected a row").with_context("cause", error.to_string())
            })?;
    }
    if group_cursor != range.group_start + range.group_count
        || atom_cursor != range.atom_start + range.atom_count
    {
        return Err(schema_error(
            "MMTF model hierarchy traversal was inconsistent",
        ));
    }
    add_global_bonds(decoded, range, &keep_map, &mut bond_builder)?;
    finish_model_data(data, builder, bond_builder, model, range)
}

fn finish_model_data(
    mut data: StructureData,
    builder: ChunkBuilder,
    bond_builder: BondTableBuilder,
    model: usize,
    range: ModelRange,
) -> Result<Structure, Diagnostic> {
    let (chunks, frame) = builder.finish();
    data.chunks = chunks.into();
    data.bonds = bond_builder.finish();
    data.topology
        .models
        .push(
            i32::try_from(model + 1).map_err(|_| schema_error("too many models"))?,
            0..u32_of(range.chain_count, "chain count")?,
        )
        .map_err(|error| {
            schema_error("model table rejected a row").with_context("cause", error.to_string())
        })?;
    data.coords = CoordinateStore::Single(frame);
    Ok(Structure::new(data))
}

fn dense_ensemble(models: &[Structure]) -> Result<Structure, Diagnostic> {
    let first = models
        .first()
        .ok_or_else(|| schema_error("MMTF contains no models"))?;
    let mut data = first.data().clone();
    data.topology.models = pdbiox_core::topology::ModelTable::default();
    let mut frames = Vec::with_capacity(models.len());
    for (position, model) in models.iter().enumerate() {
        let mut frame = CoordinateBlock::with_capacity(model.atom_count() as usize);
        for coordinate in model.positions() {
            frame.push(*coordinate);
        }
        frames.push(frame);
        data.topology
            .models
            .push(
                i32::try_from(position + 1).map_err(|_| schema_error("too many models"))?,
                0..u32_of(data.topology.chains.len(), "chain count")?,
            )
            .map_err(|error| {
                schema_error("model table rejected a row").with_context("cause", error.to_string())
            })?;
    }
    data.coords = CoordinateStore::Dense { frames };
    Ok(Structure::new(data))
}

fn ragged_ensemble(models: Vec<Structure>) -> Result<Structure, Diagnostic> {
    let mut data = StructureData::empty();
    if let Some(first) = models.first() {
        data.entry = first.data().entry.clone();
        data.cell = first.data().cell;
    }
    for position in 0..models.len() {
        let number = i32::try_from(position + 1)
            .map_err(|_| schema_error("MMTF model count exceeds signed 32-bit range"))?;
        data.topology.models.push(number, 0..0).map_err(|error| {
            schema_error("model table rejected a row").with_context("cause", error.to_string())
        })?;
    }
    data.coords = CoordinateStore::Ragged { models };
    Ok(Structure::new(data))
}

fn u32_of(value: usize, field: &'static str) -> Result<u32, Diagnostic> {
    u32::try_from(value).map_err(|_| {
        schema_error("MMTF index exceeds unsigned 32-bit range").with_context("field", field)
    })
}

include!("reader/helpers.rs");
