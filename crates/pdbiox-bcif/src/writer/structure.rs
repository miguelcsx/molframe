//! Direct structure-to-BinaryCIF orchestration.

use crate::container::EncodedCategory;
use crate::messagepack::Consuming;
use pdbiox_cif::{CanonicalProjection, CifWriteError, CifWriteOptions, canonical_projection};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES;
use pdbiox_core::structure::Structure;
use serde::Serialize;
use serde::ser::{Error as _, SerializeSeq, SerializeStruct, Serializer};
use std::cell::RefCell;
use std::io::Write;

const PLANNED_ATOM_COLUMN_BYTES_PER_ROW: usize = 40;

/// Writes a structure through the typed canonical projection to memory.
///
/// This is the explicit in-memory convenience API. File output should use
/// [`write_structure_to`] so it never retains the complete encoded file.
///
/// # Errors
///
/// Returns canonical projection or `BinaryCIF` encoding diagnostics.
pub fn write_structure(structure: &Structure) -> Result<Vec<u8>, Vec<Diagnostic>> {
    write_structure_with_options(structure, &CifWriteOptions::new())
}

/// Writes a structure with explicit canonical decisions to memory.
///
/// # Errors
///
/// Returns canonical projection or `BinaryCIF` encoding diagnostics.
pub fn write_structure_with_options(
    structure: &Structure,
    options: &CifWriteOptions,
) -> Result<Vec<u8>, Vec<Diagnostic>> {
    let mut output = Vec::new();
    write_structure_to(structure, options, &mut output)?;
    Ok(output)
}

/// Streams deterministic `BinaryCIF` to a byte destination.
///
/// Atom-site columns are generated, serialized and released one at a time.
/// Neither a row document nor a complete encoded file is retained.
///
/// # Errors
///
/// Returns canonical projection, column encoding, or destination diagnostics.
pub fn write_structure_to<W: Write>(
    structure: &Structure,
    options: &CifWriteOptions,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    write_structure_to_with_memory_limit(
        structure,
        options,
        DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES,
        output,
    )
}

/// Streams `BinaryCIF` under an explicit per-column workspace ceiling.
///
/// # Errors
///
/// Returns a memory-policy diagnostic before output when the largest planned
/// atom-column encode cannot fit, plus the normal projection/encoding errors.
pub fn write_structure_to_with_memory_limit<W: Write>(
    structure: &Structure,
    options: &CifWriteOptions,
    memory_limit_bytes: usize,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    let projection = canonical_projection(structure, options)
        .map_err(|error| vec![projection_diagnostic(&error)])?;
    let rows = projection
        .atom_site_row_count()
        .map_err(|error| vec![projection_diagnostic(&error)])?;
    validate_workspace(rows, memory_limit_bytes).map_err(|error| vec![error])?;
    let metadata = super::metadata::encode(projection).map_err(|error| vec![error])?;
    let connections = super::connections::encode(projection).map_err(|error| vec![error])?;
    let file = StreamingFile {
        projection,
        rows,
        encoder: format!("pdbiox {}", env!("CARGO_PKG_VERSION")),
        metadata: RefCell::new(Some(metadata)),
        connections: RefCell::new(connections),
    };
    rmp_serde::encode::write_named(output, &file)
        .map_err(|error| vec![Diagnostic::new(Code::E1401).with_message(error.to_string())])
}

fn validate_workspace(rows: usize, limit: usize) -> Result<(), Diagnostic> {
    let required = rows
        .checked_mul(PLANNED_ATOM_COLUMN_BYTES_PER_ROW)
        .ok_or_else(|| memory_diagnostic(usize::MAX, limit))?;
    if limit == 0 || required > limit {
        return Err(memory_diagnostic(required, limit));
    }
    Ok(())
}

fn memory_diagnostic(required: usize, limit: usize) -> Diagnostic {
    Diagnostic::new(Code::E7901)
        .with_message("BinaryCIF atom-column workspace exceeds the output policy")
        .with_context("required", required.to_string())
        .with_context("limit", limit.to_string())
}

struct StreamingFile<'a> {
    projection: CanonicalProjection<'a>,
    rows: usize,
    encoder: String,
    metadata: RefCell<Option<Vec<EncodedCategory>>>,
    connections: RefCell<Option<EncodedCategory>>,
}

impl Serialize for StreamingFile<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("EncodedFile", 3)?;
        state.serialize_field("version", "0.3.0")?;
        state.serialize_field("encoder", &self.encoder)?;
        state.serialize_field("dataBlocks", &Blocks(self))?;
        state.end()
    }
}

struct Blocks<'a, 'data>(&'a StreamingFile<'data>);

impl Serialize for Blocks<'_, '_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut sequence = serializer.serialize_seq(Some(1))?;
        sequence.serialize_element(&Block(self.0))?;
        sequence.end()
    }
}

struct Block<'a, 'data>(&'a StreamingFile<'data>);

impl Serialize for Block<'_, '_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("EncodedBlock", 2)?;
        state.serialize_field("header", self.0.projection.block_id())?;
        state.serialize_field("categories", &Categories(self.0))?;
        state.end()
    }
}

struct Categories<'a, 'data>(&'a StreamingFile<'data>);

impl Serialize for Categories<'_, '_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let metadata = self
            .0
            .metadata
            .borrow_mut()
            .take()
            .ok_or_else(|| S::Error::custom("BinaryCIF metadata was already serialized"))?;
        let connection = self.0.connections.borrow_mut().take();
        let length = metadata.len() + 1 + usize::from(connection.is_some());
        let mut sequence = serializer.serialize_seq(Some(length))?;
        for category in metadata {
            sequence.serialize_element(&Consuming::new(category))?;
        }
        sequence.serialize_element(&AtomCategory {
            projection: self.0.projection,
            rows: self.0.rows,
        })?;
        if let Some(category) = connection {
            sequence.serialize_element(&Consuming::new(category))?;
        }
        sequence.end()
    }
}

struct AtomCategory<'a> {
    projection: CanonicalProjection<'a>,
    rows: usize,
}

impl Serialize for AtomCategory<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("EncodedCategory", 3)?;
        state.serialize_field("name", "_atom_site")?;
        state.serialize_field("rowCount", &self.rows)?;
        state.serialize_field(
            "columns",
            &AtomColumns {
                projection: self.projection,
                rows: self.rows,
            },
        )?;
        state.end()
    }
}

struct AtomColumns<'a> {
    projection: CanonicalProjection<'a>,
    rows: usize,
}

impl Serialize for AtomColumns<'_> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut sequence = serializer.serialize_seq(Some(super::atoms::COLUMN_COUNT))?;
        super::atoms::visit_columns(self.projection, self.rows, |column| {
            sequence
                .serialize_element(&Consuming::new(column))
                .map_err(encoding_diagnostic)
        })
        .map_err(S::Error::custom)?;
        sequence.end()
    }
}

fn encoding_diagnostic(error: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message(error.to_string())
}

pub(super) fn projection_diagnostic(error: &CifWriteError) -> Diagnostic {
    Diagnostic::new(Code::E4105)
        .with_message("structure cannot be projected to canonical BinaryCIF")
        .with_context("reason", error.to_string())
}

#[cfg(test)]
#[path = "structure_tests.rs"]
mod tests;
