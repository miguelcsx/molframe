//! Hierarchy compaction after atom selection materialisation.

use super::StructureData;
use crate::chunk::ChunkBuilder;
use crate::diagnostic::{Code, Diagnostic};
use crate::index::{ChainIndex, EntityIndex, ModelIndex, ResidueIndex};
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::topology::{ChainRecord, ResidueRecord, Topology};

pub(super) fn compact_hierarchy(data: &mut StructureData) -> Result<(), Diagnostic> {
    let selected_residues: Vec<bool> = (0..u32_of(data.topology.residues.len())?)
        .map(ResidueIndex::new)
        .map(|residue| {
            data.topology
                .residues
                .atoms(residue)
                .is_some_and(|range| !range.is_empty())
        })
        .collect();
    let selected_chains: Vec<bool> = data
        .topology
        .chains
        .iter()
        .map(|chain| {
            data.topology.chains.residues(chain).is_some_and(|range| {
                range
                    .filter_map(|position| selected_residues.get(position as usize))
                    .any(|selected| *selected)
            })
        })
        .collect();
    let selected_entities = selected_entities(data, &selected_chains);
    let mut topology = Topology::default();
    let mut entity_remap = vec![u32::MAX; data.topology.entities.len()];
    copy_entities(data, &mut topology, &selected_entities, &mut entity_remap)?;
    let mut residue_remap = vec![u32::MAX; data.topology.residues.len()];
    copy_models_and_children(
        data,
        &mut topology,
        &selected_chains,
        &selected_residues,
        &entity_remap,
        &mut residue_remap,
    )?;
    data.chunks = rebuild_chunks(data, &residue_remap)?;
    data.topology = topology;
    Ok(())
}

fn selected_entities(data: &StructureData, selected_chains: &[bool]) -> Vec<bool> {
    let mut selected_entities = vec![false; data.topology.entities.len()];
    for chain in data.topology.chains.iter() {
        if selected_chains.get(chain.as_usize()) != Some(&true) {
            continue;
        }
        if let Some(entity) = data.topology.chains.entity(chain)
            && let Some(selected) = selected_entities.get_mut(entity.as_usize())
        {
            *selected = true;
        }
    }
    selected_entities
}

fn copy_entities(
    data: &StructureData,
    topology: &mut Topology,
    selected: &[bool],
    remap: &mut [u32],
) -> Result<(), Diagnostic> {
    for old in data.topology.entities.iter() {
        if selected.get(old.as_usize()) != Some(&true) {
            continue;
        }
        let kind = required(data.topology.entities.kind(old), "entity kind")?;
        let description = optional_symbol(data.topology.entities.description(old));
        let sequence = data.topology.entities.canonical_sequence(old);
        let new = match data.topology.entities.id(old) {
            Some(id) => topology.entities.push(id, kind, description, sequence),
            None => topology
                .entities
                .push_without_id(kind, description, sequence),
        }
        .map_err(|_| invalid_hierarchy())?;
        if let Some(slot) = remap.get_mut(old.as_usize()) {
            *slot = new.get();
        }
    }
    Ok(())
}

fn copy_models_and_children(
    data: &StructureData,
    topology: &mut Topology,
    selected_chains: &[bool],
    selected_residues: &[bool],
    entity_remap: &[u32],
    residue_remap: &mut [u32],
) -> Result<(), Diagnostic> {
    for model in data.topology.models.iter() {
        let first_chain = u32_of(topology.chains.len())?;
        let chains = required(data.topology.models.chains(model), "model chains")?;
        for chain_position in chains {
            let old_chain = ChainIndex::new(chain_position);
            if selected_chains.get(old_chain.as_usize()) != Some(&true) {
                continue;
            }
            copy_chain(
                data,
                topology,
                old_chain,
                selected_residues,
                entity_remap,
                residue_remap,
            )?;
        }
        let model_number = required(data.topology.models.model_num(model), "model number")?;
        topology
            .models
            .push(model_number, first_chain..u32_of(topology.chains.len())?)
            .map_err(|_| invalid_hierarchy())?;
    }
    Ok(())
}

fn copy_chain(
    data: &StructureData,
    topology: &mut Topology,
    old_chain: ChainIndex,
    selected_residues: &[bool],
    entity_remap: &[u32],
    residue_remap: &mut [u32],
) -> Result<(), Diagnostic> {
    let first_residue = u32_of(topology.residues.len())?;
    let residues = required(data.topology.chains.residues(old_chain), "chain residues")?;
    for residue_position in residues {
        let old_residue = ResidueIndex::new(residue_position);
        if selected_residues.get(old_residue.as_usize()) != Some(&true) {
            continue;
        }
        copy_residue(data, topology, old_residue, residue_remap)?;
    }
    let old_entity = required(data.topology.chains.entity(old_chain), "chain entity")?;
    let mapped_entity = *entity_remap
        .get(old_entity.as_usize())
        .filter(|mapped| **mapped != u32::MAX)
        .ok_or_else(invalid_hierarchy)?;
    topology
        .chains
        .push(
            ChainRecord {
                label_asym_id: required(
                    data.topology.chains.label_asym_id(old_chain),
                    "chain label",
                )?,
                auth_asym_id: optional_symbol(data.topology.chains.auth_asym_id(old_chain)),
                entity: EntityIndex::new(mapped_entity),
                polymer_kind: required(
                    data.topology.chains.polymer_kind(old_chain),
                    "polymer kind",
                )?,
            },
            first_residue..u32_of(topology.residues.len())?,
        )
        .map_err(|_| invalid_hierarchy())?;
    Ok(())
}

fn copy_residue(
    data: &StructureData,
    topology: &mut Topology,
    old: ResidueIndex,
    remap: &mut [u32],
) -> Result<(), Diagnostic> {
    let old_atoms = required(data.topology.residues.atoms(old), "residue atoms")?;
    let atom_count = old_atoms
        .end
        .checked_sub(old_atoms.start)
        .ok_or_else(invalid_hierarchy)?;
    let first_atom = topology.atom_count();
    let new = topology
        .residues
        .push(
            ResidueRecord {
                label_comp_id: required(
                    data.topology.residues.label_comp_id(old),
                    "residue component",
                )?,
                auth_comp_id: optional_symbol(data.topology.residues.auth_comp_id(old)),
                label_seq_id: optional_i32(data.topology.residues.label_seq_id(old)),
                auth_seq_id: optional_i32(data.topology.residues.auth_seq_id(old)),
                ins_code: optional_symbol(data.topology.residues.ins_code(old)),
                het: data.topology.residues.is_het(old),
            },
            first_atom..first_atom + atom_count,
        )
        .map_err(|_| invalid_hierarchy())?;
    if let Some(slot) = remap.get_mut(old.as_usize()) {
        *slot = new.get();
    }
    Ok(())
}

fn u32_of(value: usize) -> Result<u32, Diagnostic> {
    u32::try_from(value).map_err(|_| invalid_hierarchy())
}

fn rebuild_chunks(
    data: &StructureData,
    residue_remap: &[u32],
) -> Result<std::sync::Arc<Vec<crate::chunk::AtomChunk>>, Diagnostic> {
    let positions = data
        .coords
        .block(ModelIndex::new(0))
        .ok_or_else(|| Diagnostic::new(Code::E6003))?
        .as_slice();
    let mut builder = ChunkBuilder::new();
    for chunk in data.chunks.iter() {
        builder.start_model(chunk.model());
        for atom in chunk.atoms() {
            let local = atom
                .checked_sub(chunk.atoms().start)
                .ok_or_else(invalid_hierarchy)?;
            let mut record = chunk
                .record(
                    local,
                    &data.topology.residues,
                    positions.get(atom as usize).copied(),
                )
                .ok_or_else(invalid_hierarchy)?;
            let mapped = *residue_remap
                .get(record.residue.as_usize())
                .filter(|value| **value != u32::MAX)
                .ok_or_else(invalid_hierarchy)?;
            record.residue = ResidueIndex::new(mapped);
            builder.push(record);
        }
    }
    let (chunks, _) = builder.finish();
    Ok(std::sync::Arc::new(chunks))
}

fn optional_symbol(value: Option<crate::SymbolId>) -> OptionalSymbol {
    value.map_or(OptionalSymbol::NONE, OptionalSymbol::some)
}

fn optional_i32(value: Option<i32>) -> OptionalI32 {
    value.map_or(OptionalI32::NONE, OptionalI32::some)
}

fn required<T>(value: Option<T>, field: &'static str) -> Result<T, Diagnostic> {
    value.ok_or_else(|| invalid_hierarchy().with_context("field", field))
}

fn invalid_hierarchy() -> Diagnostic {
    Diagnostic::new(Code::E3001)
}
