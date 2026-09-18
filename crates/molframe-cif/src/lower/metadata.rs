//! Reconcile entry and entity metadata collected after coordinates started.

use super::entry::{AsymEntity, read_asym_entities, read_cell, read_entities, read_entry};
use crate::document::DataBlock;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::index::EntityIndex;
use molframe_core::io::ReadOptions;
use molframe_core::optional::OptionalSymbol;
use molframe_core::structure::StructureData;
use molframe_core::symbol::SymbolId;
use molframe_core::topology::{ChainRecord, ChainTable, EntityKind, EntityTable, PolymerKind};
use std::ops::Range;

pub(super) fn refresh(
    block: &DataBlock,
    options: &ReadOptions,
    data: &mut StructureData,
    findings: &mut Diagnostics,
) {
    if !options.only_atomic_coords {
        data.entry = molframe_core::structure::EntryMetadata::default();
        data.cell = None;
        read_entry(block, data, findings);
        read_cell(block, data);
    }
    refresh_entities(block, data, findings);
}

#[derive(Clone)]
struct ChainSnapshot {
    label: SymbolId,
    auth: OptionalSymbol,
    entity_id: Option<SymbolId>,
    residues: Range<u32>,
}

fn refresh_entities(block: &DataBlock, data: &mut StructureData, findings: &mut Diagnostics) {
    let chains = chain_snapshots(data);
    data.topology.entities = EntityTable::default();
    read_entities(block, data, findings);
    let mappings = read_asym_entities(block, data, findings);
    retain_provisional_entities(&chains, data, findings);

    let mut rebuilt = ChainTable::default();
    for chain in chains {
        let Some(entity) = resolve_entity(&chain, &mappings, data, findings) else {
            continue;
        };
        let polymer_kind = match data.topology.entities.kind(entity) {
            Some(EntityKind::Polymer) => PolymerKind::Other,
            _ => PolymerKind::None,
        };
        if rebuilt
            .push(
                ChainRecord {
                    label_asym_id: chain.label,
                    auth_asym_id: chain.auth,
                    entity,
                    polymer_kind,
                },
                chain.residues,
            )
            .is_err()
        {
            findings.push(Diagnostic::new(Code::E3001));
        }
    }
    data.topology.chains = rebuilt;
}

fn chain_snapshots(data: &StructureData) -> Vec<ChainSnapshot> {
    let mut snapshots = Vec::with_capacity(data.topology.chains.len());
    for chain in data.topology.chains.iter() {
        let (Some(label), Some(residues), Some(entity)) = (
            data.topology.chains.label_asym_id(chain),
            data.topology.chains.residues(chain),
            data.topology.chains.entity(chain),
        ) else {
            continue;
        };
        snapshots.push(ChainSnapshot {
            label,
            auth: OptionalSymbol::from(data.topology.chains.auth_asym_id(chain)),
            entity_id: data.topology.entities.id(entity),
            residues,
        });
    }
    snapshots
}

fn retain_provisional_entities(
    chains: &[ChainSnapshot],
    data: &mut StructureData,
    findings: &mut Diagnostics,
) {
    for id in chains.iter().filter_map(|chain| chain.entity_id) {
        if data.topology.entities.find_by_id(id).is_some() {
            continue;
        }
        if data
            .topology
            .entities
            .push(id, EntityKind::Unknown, OptionalSymbol::NONE, &[])
            .is_err()
        {
            findings.push(Diagnostic::new(Code::E3001));
        }
    }
}

fn resolve_entity(
    chain: &ChainSnapshot,
    mappings: &[AsymEntity],
    data: &mut StructureData,
    findings: &mut Diagnostics,
) -> Option<EntityIndex> {
    if let Some(id) = chain.entity_id
        && let Some(entity) = data.topology.entities.find_by_id(id)
    {
        return Some(entity);
    }
    if let Some(mapping) = mappings
        .iter()
        .find(|mapping| mapping.asym_id == chain.label)
    {
        return Some(mapping.entity);
    }
    if let Some(entity) = data
        .topology
        .entities
        .iter()
        .find(|entity| data.topology.entities.id(*entity).is_none())
    {
        return Some(entity);
    }
    match data
        .topology
        .entities
        .push_without_id(EntityKind::Unknown, OptionalSymbol::NONE, &[])
    {
        Ok(entity) => Some(entity),
        Err(error) => {
            findings.push(Diagnostic::new(Code::E3001).with_context("cause", error.to_string()));
            None
        }
    }
}
