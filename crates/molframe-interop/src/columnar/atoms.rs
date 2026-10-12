//! Arrow atom-table batches aligned to molframe chunk boundaries.

use super::dictionary::{AtomDictionaries, dictionary_type, present};
use super::extension::{ExportCost, field};
use super::owner::{SnapshotOwner, f32_buffer};
use super::stream::{ArrowStream, ArrowTableExport};
use crate::numeric::f64_to_f32;
use arrow::array::{
    ArrayRef, BooleanArray, BooleanBuilder, FixedSizeListArray, Float32Array, UInt8Array,
    UInt32Array, builder::NullBufferBuilder,
};
use arrow::buffer::NullBuffer;
use arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use arrow::error::{ArrowError, Result};
use arrow::record_batch::RecordBatch;
use molframe_core::annotation::{AtomAnnotation, PLDDT_ANNOTATION};
use molframe_core::topology::ResidueTable;
use molframe_core::{AtomChunk, ResidueIndex, Structure};
use std::collections::HashMap;
use std::ops::Range;
use std::sync::{Arc, OnceLock};

/// Schema metadata key carrying the atom-table schema version.
pub(crate) const SCHEMA_VERSION_KEY: &str = "molframe:schema_version";
/// Current atom-table schema version.
pub(crate) const SCHEMA_VERSION: &str = "2";
const UNIT_KEY: &str = "molframe:unit";

/// Arrow view of a structure's atom table.
///
/// Cloning is `O(1)`: the immutable snapshot is shared through [`Arc`] and the
/// zero-copy buffers alias it, so a clone bumps reference counts rather than
/// copying coordinates.
#[derive(Clone, Debug)]
pub struct AtomTable {
    structure: Arc<Structure>,
    schema: SchemaRef,
    owner: SnapshotOwner,
    dictionaries: Arc<OnceLock<AtomDictionaries>>,
}

impl AtomTable {
    /// Binds the immutable structure snapshot to an Arrow table.
    #[must_use]
    pub fn new(structure: &Structure) -> Self {
        let structure = Arc::new(structure.clone());
        Self {
            structure: structure.clone(),
            schema: Arc::new(atom_schema(&structure)),
            owner: SnapshotOwner::from_arc(structure),
            dictionaries: Arc::new(OnceLock::new()),
        }
    }

    /// Categorical dictionaries, built on first use and shared by every batch.
    fn dictionaries(&self) -> &AtomDictionaries {
        self.dictionaries
            .get_or_init(|| AtomDictionaries::build(&self.structure))
    }

    /// Arrow schema with molframe extension and copy-cost metadata.
    #[must_use]
    pub fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    /// One record batch per internal atom chunk.
    ///
    /// # Errors
    ///
    /// Returns an Arrow error if a chunk violates the declared schema or a
    /// backing allocation cannot be represented.
    pub fn record_batches(&self) -> Result<Vec<RecordBatch>> {
        <Self as ArrowTableExport>::record_batches(self)
    }

    /// Creates a consumable Arrow C Stream Interface value.
    ///
    /// Construction is constant-space: chunks are converted only when the
    /// consumer pulls them from the stream. A conversion error is reported by
    /// that pull rather than by this method.
    ///
    /// # Errors
    ///
    /// Reserved for failures that occur while constructing the stream adapter.
    pub fn arrow_stream(&self) -> Result<ArrowStream> {
        <Self as ArrowTableExport>::arrow_stream(self)
    }

    /// As [`Self::arrow_stream`], restricted to the atoms in `rows`.
    ///
    /// Chunks wholly outside the window are never converted.
    ///
    /// # Errors
    ///
    /// Returns an Arrow error when the window lies outside the table.
    pub fn arrow_stream_rows(&self, rows: std::ops::Range<usize>) -> Result<ArrowStream> {
        <Self as ArrowTableExport>::arrow_stream_rows(self, rows)
    }

    fn chunk_batch(&self, chunk: &AtomChunk) -> Result<RecordBatch> {
        let range = chunk.atoms();
        let positions = self
            .structure
            .positions()
            .get(range.start as usize..range.end as usize)
            .ok_or_else(|| {
                ArrowError::InvalidArgumentError("chunk coordinates absent".to_owned())
            })?;
        let occupancies = chunk.occupancies_plain().ok_or_else(|| {
            ArrowError::InvalidArgumentError("occupancies require decode".to_owned())
        })?;
        let b_factors = chunk.b_factors_plain().ok_or_else(|| {
            ArrowError::InvalidArgumentError("temperature factors require decode".to_owned())
        })?;
        let index = UInt32Array::from_iter_values(range.clone());
        let residue = residue_indices(range.clone(), &self.structure.data().topology.residues)?;
        let element = UInt8Array::from_iter_values((0..chunk.len()).filter_map(|local| {
            chunk
                .element(local)
                .map(molframe_core::Element::atomic_number)
        }));
        let coordinates = coordinates_array(positions, chunk, &self.owner)?;
        let dictionaries = self.dictionaries();
        let residue_key = |keys: &[i32]| {
            residue
                .values()
                .iter()
                .map(|raw| keys.get(*raw as usize).copied().and_then(present))
                .collect::<Vec<_>>()
        };
        let chain = dictionaries
            .chain
            .array(residue_key(&dictionaries.residue_chain))?;
        let label_asym = dictionaries
            .label_asym
            .array(residue_key(&dictionaries.residue_label_asym))?;
        let resname = dictionaries
            .resname
            .array(residue_key(&dictionaries.residue_resname))?;
        let names = dictionaries
            .atom_name
            .array(match chunk.atom_names_plain() {
                Some(names) => names
                    .iter()
                    .map(|symbol| dictionaries.atom_name.key(symbol.get()))
                    .collect::<Vec<_>>(),
                None => (0..chunk.len())
                    .map(|local| {
                        chunk
                            .atom_name(local)
                            .and_then(|symbol| dictionaries.atom_name.key(symbol.get()))
                    })
                    .collect(),
            })?;
        let element_symbol = dictionaries.element.array((0..chunk.len()).map(|local| {
            chunk
                .element(local)
                .and_then(|element| dictionaries.element.key(u32::from(element.atomic_number())))
        }))?;
        let altloc = dictionaries.altloc.array((0..chunk.len()).map(|local| {
            chunk
                .alt_id(local)
                .filter(|alt| !alt.is_blank())
                .and_then(|alt| dictionaries.altloc.key(alt.get()))
        }))?;
        let (occupancy_validity, occupancy_unknown) =
            presence_arrays(chunk.len(), |local| chunk.occupancy(local));
        let occupancy =
            Float32Array::new(f32_buffer(occupancies, &self.owner)?, occupancy_validity);
        let (b_factor_validity, b_factor_unknown) =
            presence_arrays(chunk.len(), |local| chunk.b_factor(local));
        let b_factor = Float32Array::new(f32_buffer(b_factors, &self.owner)?, b_factor_validity);
        let mut arrays: Vec<ArrayRef> = vec![
            Arc::new(index),
            Arc::new(residue),
            chain,
            label_asym,
            resname,
            Arc::new(coordinates),
            Arc::new(element),
            element_symbol,
            names,
            altloc,
            Arc::new(occupancy),
            Arc::new(occupancy_unknown),
            Arc::new(b_factor),
            Arc::new(b_factor_unknown),
        ];
        if let Some(column) = plddt_column(&self.structure) {
            arrays.push(Arc::new(Float32Array::from_iter(
                (range.start..range.end).map(|atom| {
                    column
                        .get(atom)
                        .and_then(|(value, state)| state.is_present().then(|| f64_to_f32(value)))
                }),
            )));
        }
        RecordBatch::try_new(self.schema.clone(), arrays)
    }
}

/// Expands ordered residue ranges in `O(atoms + overlapping residues)` time.
fn residue_indices(range: Range<u32>, residues: &ResidueTable) -> Result<UInt32Array> {
    let atom_count = range.end.checked_sub(range.start).ok_or_else(|| {
        ArrowError::InvalidArgumentError("atom chunk range is reversed".to_owned())
    })?;
    let expected = usize::try_from(atom_count)
        .map_err(|_| ArrowError::MemoryError("atom chunk length exceeds usize".to_owned()))?;
    if expected == 0 {
        return Ok(UInt32Array::from(Vec::<u32>::new()));
    }
    let last_atom = range
        .end
        .checked_sub(1)
        .ok_or_else(|| ArrowError::InvalidArgumentError("empty atom chunk".to_owned()))?;
    let first_residue = residues.containing(range.start).ok_or_else(|| {
        ArrowError::InvalidArgumentError("first atom has no containing residue".to_owned())
    })?;
    let last_residue = residues.containing(last_atom).ok_or_else(|| {
        ArrowError::InvalidArgumentError("last atom has no containing residue".to_owned())
    })?;
    let mut values = Vec::with_capacity(expected);
    for raw in first_residue.get()..=last_residue.get() {
        let atoms = residues.atoms(ResidueIndex::new(raw)).ok_or_else(|| {
            ArrowError::InvalidArgumentError("residue atom range absent".to_owned())
        })?;
        let overlap_start = atoms.start.max(range.start);
        let overlap_end = atoms.end.min(range.end);
        let overlap = usize::try_from(overlap_end.saturating_sub(overlap_start))
            .map_err(|_| ArrowError::MemoryError("residue overlap exceeds usize".to_owned()))?;
        let next = values
            .len()
            .checked_add(overlap)
            .ok_or_else(|| ArrowError::MemoryError("residue index length overflow".to_owned()))?;
        if next > expected {
            return Err(ArrowError::InvalidArgumentError(
                "residue ranges overlap within atom chunk".to_owned(),
            ));
        }
        values.resize(next, raw);
    }
    if values.len() != expected {
        return Err(ArrowError::InvalidArgumentError(
            "residue ranges do not cover atom chunk".to_owned(),
        ));
    }
    Ok(UInt32Array::from(values))
}

impl ArrowTableExport for AtomTable {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn batch_count(&self) -> usize {
        self.structure.data().chunks.len()
    }

    fn row_count(&self) -> usize {
        self.structure.atom_count() as usize
    }

    fn batch_rows(&self, index: usize) -> std::ops::Range<usize> {
        match self.structure.data().chunks.get(index) {
            Some(chunk) => chunk.atoms().start as usize..chunk.atoms().end as usize,
            None => 0..0,
        }
    }

    fn batch(&self, index: usize) -> Result<RecordBatch> {
        let chunk = self.structure.data().chunks.get(index).ok_or_else(|| {
            ArrowError::InvalidArgumentError("atom chunk index is out of bounds".to_owned())
        })?;
        self.chunk_batch(chunk)
    }
}

fn coordinates_array(
    positions: &[[f32; 3]],
    chunk: &AtomChunk,
    owner: &SnapshotOwner,
) -> Result<FixedSizeListArray> {
    let count = positions
        .len()
        .checked_mul(3)
        .ok_or_else(|| ArrowError::MemoryError("coordinate length overflow".to_owned()))?;
    // SAFETY: `[f32; 3]` is three contiguous, equally aligned `f32` values with
    // no padding. `count` is the checked scalar length of the live input slice.
    let values = unsafe { std::slice::from_raw_parts(positions.as_ptr().cast::<f32>(), count) };
    let values: ArrayRef = Arc::new(Float32Array::new(f32_buffer(values, owner)?, None));
    FixedSizeListArray::try_new(
        Arc::new(Field::new("item", DataType::Float32, false)),
        3,
        values,
        validity(chunk.len(), |local| {
            chunk
                .has_position(local)
                .then_some(((), molframe_core::Presence::Present))
        }),
    )
}

fn validity<T>(
    len: u32,
    mut value: impl FnMut(u32) -> Option<(T, molframe_core::Presence)>,
) -> Option<NullBuffer> {
    let mut validity = NullBufferBuilder::new(len as usize);
    for position in 0..len {
        let present = value(position).is_some_and(|(_, state)| state.is_present());
        validity.append(present);
    }
    validity.finish()
}

fn presence_arrays<T>(
    len: u32,
    mut value: impl FnMut(u32) -> Option<(T, molframe_core::Presence)>,
) -> (Option<NullBuffer>, BooleanArray) {
    let mut validity = NullBufferBuilder::new(len as usize);
    let mut unknown = BooleanBuilder::with_capacity(len as usize);
    for position in 0..len {
        let state = value(position).map(|(_, state)| state);
        let present = state.is_some_and(molframe_core::Presence::is_present);
        validity.append(present);
        unknown.append_value(state == Some(molframe_core::Presence::Unknown));
    }
    (validity.finish(), unknown.finish())
}

/// The lowered pLDDT confidence column, when the structure carries one.
fn plddt_column(
    structure: &Structure,
) -> Option<&molframe_core::annotation::AnnotationColumn<f64>> {
    match structure.annotations().get(PLDDT_ANNOTATION) {
        Some(AtomAnnotation::Real(column)) => Some(column),
        _ => None,
    }
}

/// Schema-level provenance, units and version; merged under caller metadata.
fn schema_metadata(structure: &Structure) -> HashMap<String, String> {
    let entry = &structure.data().entry;
    let mut metadata = HashMap::from([
        (SCHEMA_VERSION_KEY.to_owned(), SCHEMA_VERSION.to_owned()),
        (
            "molframe:producer".to_owned(),
            format!("molframe-interop {}", env!("CARGO_PKG_VERSION")),
        ),
        ("molframe:table".to_owned(), "atoms".to_owned()),
        (
            "molframe:atom_count".to_owned(),
            structure.atom_count().to_string(),
        ),
        (
            "molframe:units".to_owned(),
            concat!(
                r#"{"coordinates":"angstrom","b_factor":"angstrom^2","#,
                r#""occupancy":"fraction","plddt":"score_0_100"}"#
            )
            .to_owned(),
        ),
    ]);
    let optional = [
        ("molframe:entry_id", entry.id.as_deref()),
        ("molframe:method", entry.method.as_deref()),
        ("molframe:input_name", entry.input_name.as_deref()),
    ];
    for (key, value) in optional {
        if let Some(value) = value {
            metadata.insert(key.to_owned(), value.to_owned());
        }
    }
    if let Some(digest) = &entry.input_sha256 {
        metadata.insert("molframe:input_sha256".to_owned(), digest.to_string());
    }
    metadata
}

fn unit(field: Field, unit: &str) -> Field {
    let mut metadata = field.metadata().clone();
    metadata.insert(UNIT_KEY.to_owned(), unit.to_owned());
    field.with_metadata(metadata)
}

fn atom_schema(structure: &Structure) -> Schema {
    let coordinates =
        DataType::FixedSizeList(Arc::new(Field::new("item", DataType::Float32, false)), 3);
    let categorical = |name: &str| field(name, dictionary_type(), true, None, ExportCost::Decode);
    let mut fields = vec![
        field(
            "atom_index",
            DataType::UInt32,
            false,
            Some("molframe.atom_index"),
            ExportCost::Decode,
        ),
        field(
            "residue_index",
            DataType::UInt32,
            false,
            Some("molframe.residue_index"),
            ExportCost::Decode,
        ),
        categorical("chain_id"),
        categorical("label_asym_id"),
        categorical("resname"),
        unit(
            field(
                "coordinates",
                coordinates,
                true,
                Some("molframe.coordinates3f"),
                ExportCost::ZeroCopy,
            ),
            "angstrom",
        ),
        field(
            "element",
            DataType::UInt8,
            false,
            Some("molframe.element"),
            ExportCost::Decode,
        ),
        categorical("element_symbol"),
        categorical("atom_name"),
        categorical("altloc"),
        unit(
            field(
                "occupancy",
                DataType::Float32,
                true,
                None,
                ExportCost::ZeroCopy,
            ),
            "fraction",
        ),
        field(
            "occupancy_unknown",
            DataType::Boolean,
            false,
            Some("molframe.validity"),
            ExportCost::Decode,
        ),
        unit(
            field(
                "b_factor",
                DataType::Float32,
                true,
                None,
                ExportCost::ZeroCopy,
            ),
            "angstrom^2",
        ),
        field(
            "b_factor_unknown",
            DataType::Boolean,
            false,
            Some("molframe.validity"),
            ExportCost::Decode,
        ),
    ];
    if plddt_column(structure).is_some() {
        fields.push(unit(
            field("plddt", DataType::Float32, true, None, ExportCost::Decode),
            "score_0_100",
        ));
    }
    Schema::new_with_metadata(fields, schema_metadata(structure))
}

#[cfg(test)]
#[path = "atoms_tests.rs"]
mod tests;
