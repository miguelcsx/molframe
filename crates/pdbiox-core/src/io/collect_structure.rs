//! Explicit budgeted materialisation from structural batches.

use crate::chunk::{AtomRecord, ChunkBuilder};
use crate::execution::{Backpressure, BatchDemand, BatchSource, MemoryReservation};
use crate::index::{EntityIndex, ResidueIndex};
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::structure::{CoordinateStore, Structure, StructureData};
use crate::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};
use crate::{
    AltId, Batch, Code, Diagnostic, Diagnostics, Element, ExecutionContext, StructureBatch,
    StructureBatchError, SymbolId,
};

const COLLECTOR_BYTES_PER_ROW: usize = 192;
const COLLECTOR_BASE_BYTES: usize = 64 * 1024;
const PULL_BYTES: usize = 8 * 1024 * 1024;
const PULL_ROWS: usize = 65_536;

/// Materialises a bounded snapshot from a structural batch source.
///
/// Multi-model inputs become explicit ragged model snapshots. This avoids
/// retaining both a topology template and all later model rows while the source
/// is still active.
///
/// # Errors
///
/// Returns a typed source, addressability or budget error before growing the
/// collector beyond the shared execution budget.
pub fn collect_structure<S>(
    source: &mut S,
    context: &ExecutionContext,
) -> Result<(Structure, Vec<Diagnostic>), StructureBatchError>
where
    S: BatchSource<Batch = StructureBatch, Error = StructureBatchError>,
{
    let mut collector = Collector::new(context)?;
    loop {
        let available = context
            .memory_budget()
            .bytes()
            .saturating_sub(context.reserved_bytes());
        if available == 0 {
            return Err(exhausted(context, COLLECTOR_BYTES_PER_ROW));
        }
        let demand = BatchDemand::new(PULL_ROWS, PULL_BYTES.min(available));
        match source.next_batch(demand, context)? {
            Backpressure::Ready(lease) => collector.append(lease.batch())?,
            Backpressure::Pending => return Err(exhausted(context, COLLECTOR_BYTES_PER_ROW)),
            Backpressure::Finished => break,
        }
    }
    collector.finish()
}

struct Collector {
    models: Vec<Structure>,
    current: Option<ModelCollector>,
    diagnostics: Diagnostics,
    reservation: MemoryReservation,
}

impl Collector {
    fn new(context: &ExecutionContext) -> Result<Self, StructureBatchError> {
        Ok(Self {
            models: Vec::new(),
            current: None,
            diagnostics: Diagnostics::new(),
            reservation: context.try_reserve(COLLECTOR_BASE_BYTES)?,
        })
    }

    fn append(&mut self, batch: &StructureBatch) -> Result<(), StructureBatchError> {
        let growth = batch
            .rows()
            .checked_mul(COLLECTOR_BYTES_PER_ROW)
            .ok_or_else(identity_overflow)?;
        self.reservation.try_grow(growth)?;
        self.diagnostics.extend(batch.diagnostics().iter().cloned());
        for row in 0..batch.rows() {
            let model = batch.models()[row];
            if self
                .current
                .as_ref()
                .is_some_and(|current| current.model != model)
            {
                self.close_model()?;
            }
            if self.current.is_none() {
                self.current = Some(ModelCollector::new(model));
            }
            let Some(current) = &mut self.current else {
                return Err(identity_overflow());
            };
            current.push(batch, row)?;
        }
        Ok(())
    }

    fn close_model(&mut self) -> Result<(), StructureBatchError> {
        let Some(current) = self.current.take() else {
            return Ok(());
        };
        self.models.push(current.finish()?);
        Ok(())
    }

    fn finish(mut self) -> Result<(Structure, Vec<Diagnostic>), StructureBatchError> {
        self.close_model()?;
        let structure = match self.models.len() {
            0 => Structure::new(StructureData::empty()),
            1 => match self.models.pop() {
                Some(model) => model,
                None => Structure::new(StructureData::empty()),
            },
            _ => ragged_structure(self.models)?,
        };
        Ok((
            structure.retain_reservation(std::sync::Arc::new(self.reservation)),
            self.diagnostics.finish(),
        ))
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct ResidueKey {
    chain: SymbolId,
    component: SymbolId,
    sequence: i32,
    insertion: OptionalSymbol,
}

struct ModelCollector {
    model: i32,
    data: StructureData,
    chunks: ChunkBuilder,
    entity: Option<EntityIndex>,
    chain: Option<SymbolId>,
    chain_first_residue: u32,
    residue: Option<ResidueKey>,
    residue_position: u32,
    atom_position: u32,
}

impl ModelCollector {
    fn new(model: i32) -> Self {
        Self {
            model,
            data: StructureData::empty(),
            chunks: ChunkBuilder::new(),
            entity: None,
            chain: None,
            chain_first_residue: 0,
            residue: None,
            residue_position: 0,
            atom_position: 0,
        }
    }

    fn push(&mut self, batch: &StructureBatch, row: usize) -> Result<(), StructureBatchError> {
        let chain = remap(batch, batch.chains()[row], &mut self.data)?;
        let component = remap(batch, batch.components()[row], &mut self.data)?;
        let insertion = remap_optional(batch, batch.insertions()[row], &mut self.data)?;
        let key = ResidueKey {
            chain,
            component,
            sequence: batch.sequences()[row],
            insertion,
        };
        self.begin_hierarchy(key, batch.heterogens()[row] != 0)?;
        let residue = self
            .residue_position
            .checked_sub(1)
            .map(ResidueIndex::new)
            .ok_or_else(identity_overflow)?;
        let atom_name = remap(batch, batch.atoms()[row], &mut self.data)?;
        let alternate = remap_alt(batch, batch.alternates()[row], &mut self.data)?;
        let (occupancy, occupancy_presence) = batch.occupancies();
        let (b_factor, b_factor_presence) = batch.b_factors();
        let (formal_charge, formal_charge_presence) = batch.formal_charges();
        let position = if batch.coordinate_presence()[row].is_present() {
            Some(batch.positions()[row])
        } else {
            None
        };
        self.chunks.push(AtomRecord {
            position,
            element: Element::from_atomic_number(batch.elements()[row]),
            atom_name,
            auth_atom_name: OptionalSymbol::some(atom_name),
            alternate_component_id: OptionalSymbol::NONE,
            alt_id: alternate,
            residue,
            occupancy: (occupancy[row], occupancy_presence[row]),
            b_factor: (b_factor[row], b_factor_presence[row]),
            formal_charge: (formal_charge[row], formal_charge_presence[row]),
            atom_site_id: batch.atom_site_ids()[row],
        });
        self.atom_position = self
            .atom_position
            .checked_add(1)
            .ok_or_else(identity_overflow)?;
        Ok(())
    }

    fn begin_hierarchy(
        &mut self,
        key: ResidueKey,
        heterogen: bool,
    ) -> Result<(), StructureBatchError> {
        if self.chain != Some(key.chain) {
            self.close_chain()?;
            self.chain = Some(key.chain);
            self.chain_first_residue = self.residue_position;
        }
        if self.residue == Some(key) {
            return Ok(());
        }
        self.close_residue()?;
        self.data
            .topology
            .residues
            .push(
                ResidueRecord {
                    label_comp_id: key.component,
                    auth_comp_id: OptionalSymbol::some(key.component),
                    label_seq_id: optional_sequence(key.sequence),
                    auth_seq_id: optional_sequence(key.sequence),
                    ins_code: key.insertion,
                    het: heterogen,
                },
                self.atom_position..self.atom_position,
            )
            .map_err(table_error)?;
        self.residue = Some(key);
        self.residue_position = self
            .residue_position
            .checked_add(1)
            .ok_or_else(identity_overflow)?;
        Ok(())
    }

    fn close_residue(&mut self) -> Result<(), StructureBatchError> {
        let Some(position) = self.residue_position.checked_sub(1) else {
            return Ok(());
        };
        let residue = ResidueIndex::new(position);
        let start = self
            .data
            .topology
            .residues
            .atoms(residue)
            .map(|range| range.start)
            .ok_or_else(identity_overflow)?;
        self.data
            .topology
            .residues
            .set_atoms(residue, start..self.atom_position)
            .map_err(table_error)
    }

    fn close_chain(&mut self) -> Result<(), StructureBatchError> {
        self.close_residue()?;
        let Some(chain) = self.chain.take() else {
            return Ok(());
        };
        if self.chain_first_residue == self.residue_position {
            return Ok(());
        }
        let entity = self.entity()?;
        self.data
            .topology
            .chains
            .push(
                ChainRecord {
                    label_asym_id: chain,
                    auth_asym_id: OptionalSymbol::some(chain),
                    entity,
                    polymer_kind: PolymerKind::None,
                },
                self.chain_first_residue..self.residue_position,
            )
            .map_err(table_error)?;
        self.residue = None;
        Ok(())
    }

    fn entity(&mut self) -> Result<EntityIndex, StructureBatchError> {
        if let Some(entity) = self.entity {
            return Ok(entity);
        }
        let identifier = self
            .data
            .dictionary
            .intern("1")
            .map_err(|_| StructureBatchError::DictionaryFull)?;
        let entity = self
            .data
            .topology
            .entities
            .push(identifier, EntityKind::Unknown, OptionalSymbol::NONE, &[])
            .map_err(table_error)?;
        self.entity = Some(entity);
        Ok(entity)
    }

    fn finish(mut self) -> Result<Structure, StructureBatchError> {
        self.close_chain()?;
        let chain_count =
            u32::try_from(self.data.topology.chains.len()).map_err(|_| identity_overflow())?;
        self.data
            .topology
            .models
            .push(self.model, 0..chain_count)
            .map_err(table_error)?;
        let (chunks, coordinates) = self.chunks.finish();
        self.data.chunks = chunks.into();
        self.data.coords = CoordinateStore::Single(coordinates);
        Ok(Structure::new(self.data))
    }
}

fn remap(
    batch: &StructureBatch,
    symbol: SymbolId,
    data: &mut StructureData,
) -> Result<SymbolId, StructureBatchError> {
    let text = batch
        .dictionary()
        .resolve(symbol)
        .ok_or_else(identity_overflow)?;
    data.dictionary
        .intern(text)
        .map_err(|_| StructureBatchError::DictionaryFull)
}

fn remap_optional(
    batch: &StructureBatch,
    symbol: OptionalSymbol,
    data: &mut StructureData,
) -> Result<OptionalSymbol, StructureBatchError> {
    match symbol.get() {
        Some(symbol) => remap(batch, symbol, data).map(OptionalSymbol::some),
        None => Ok(OptionalSymbol::NONE),
    }
}

fn remap_alt(
    batch: &StructureBatch,
    alternate: AltId,
    data: &mut StructureData,
) -> Result<AltId, StructureBatchError> {
    let Some(symbol) = alternate.symbol() else {
        return Ok(AltId::BLANK);
    };
    AltId::labelled(remap(batch, symbol, data)?).ok_or_else(identity_overflow)
}

fn optional_sequence(value: i32) -> OptionalI32 {
    if value == i32::MIN {
        OptionalI32::NONE
    } else {
        OptionalI32::some(value)
    }
}

fn ragged_structure(models: Vec<Structure>) -> Result<Structure, StructureBatchError> {
    let mut data = StructureData::empty();
    for model in &models {
        let number = model
            .data()
            .topology
            .models
            .model_num(crate::ModelIndex::new(0))
            .ok_or_else(identity_overflow)?;
        data.topology
            .models
            .push(number, 0..0)
            .map_err(table_error)?;
    }
    data.coords = CoordinateStore::Ragged { models };
    Ok(Structure::new(data))
}

fn exhausted(context: &ExecutionContext, required: usize) -> StructureBatchError {
    StructureBatchError::RecordExceedsBudget {
        required,
        available: context
            .memory_budget()
            .bytes()
            .saturating_sub(context.reserved_bytes()),
    }
}

fn table_error(error: impl std::fmt::Display) -> StructureBatchError {
    StructureBatchError::Diagnostic(
        Diagnostic::new(Code::E3001).with_context("cause", error.to_string()),
    )
}

fn identity_overflow() -> StructureBatchError {
    StructureBatchError::Diagnostic(
        Diagnostic::new(Code::E1903).with_message("collected structure exceeds local index space"),
    )
}

#[cfg(test)]
#[path = "collect_structure_tests.rs"]
mod tests;
