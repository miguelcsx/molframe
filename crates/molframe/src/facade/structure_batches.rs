//! Public bounded structural-reader dispatch.

use molframe_core::{
    Backpressure, BatchDemand, BatchLease, BatchSource, ChunkId, Code, DatasetId, Diagnostic,
    ExecutionContext, Format, LogicalRow, ReadOptions, SourceBytes, StructureBatch,
    StructureBatchError, WindowedSourceFile,
};
use std::path::Path;

const DEFAULT_SOURCE_WINDOW_BYTES: usize = 4 * 1024 * 1024;
const MINIMUM_SOURCE_WINDOW_BYTES: usize = 4 * 1024;

/// One linked structural batch reader selected by content or explicit format.
#[derive(Debug)]
pub enum StructureBatchReader {
    /// `BinaryCIF` source.
    #[cfg(feature = "bcif")]
    BinaryCif(Box<molframe_bcif::BcifBatchSource<WindowedSourceFile>>),
    /// Text PDB source.
    #[cfg(feature = "pdb")]
    Pdb(Box<molframe_pdb::PdbBatchSource<WindowedSourceFile>>),
    /// Text mmCIF source.
    #[cfg(feature = "mmcif")]
    Mmcif(Box<molframe_cif::MmcifBatchSource<WindowedSourceFile>>),
    /// Text `ModelCIF` source selected explicitly by its constructor.
    #[cfg(feature = "modelcif")]
    ModelCif(Box<molframe_modelcif::ModelCifBatchSource<WindowedSourceFile>>),
}

impl BatchSource for StructureBatchReader {
    type Batch = StructureBatch;
    type Error = StructureBatchError;

    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<Self::Batch>>, Self::Error> {
        match self {
            #[cfg(feature = "bcif")]
            Self::BinaryCif(source) => source.next_batch(demand, context),
            #[cfg(feature = "pdb")]
            Self::Pdb(source) => source.next_batch(demand, context),
            #[cfg(feature = "mmcif")]
            Self::Mmcif(source) => source.next_batch(demand, context),
            #[cfg(feature = "modelcif")]
            Self::ModelCif(source) => source.next_batch(demand, context),
        }
    }
}

/// Opens a structural file without reserving memory proportional to file size.
///
/// # Errors
///
/// Returns a typed format, source, addressability or budget error. Formats that
/// do not yet have a linked bounded reader are refused rather than routed
/// through a whole-input compatibility path.
pub fn open_structure_batches(
    path: impl AsRef<Path>,
    options: &ReadOptions,
    context: &ExecutionContext,
) -> Result<StructureBatchReader, StructureBatchError> {
    let path = path.as_ref();
    let window_bytes = source_window_bytes(context)?;
    let mut source = WindowedSourceFile::open(path, window_bytes, context)?;
    let format = detect_format(&mut source, path, options.format, window_bytes)?;
    match format {
        #[cfg(feature = "pdb")]
        Format::Pdb => Ok(StructureBatchReader::Pdb(Box::new(
            molframe_pdb::PdbBatchSource::new(
                source,
                options.clone(),
                DatasetId::new(0),
                ChunkId::new(0),
                LogicalRow::new(0),
                window_bytes,
            )?,
        ))),
        #[cfg(feature = "mmcif")]
        Format::Mmcif => Ok(StructureBatchReader::Mmcif(Box::new(
            molframe_cif::MmcifBatchSource::new(
                source,
                options.clone(),
                DatasetId::new(0),
                ChunkId::new(0),
                LogicalRow::new(0),
                window_bytes,
                context,
            )?,
        ))),
        #[cfg(feature = "bcif")]
        Format::BinaryCif => Ok(StructureBatchReader::BinaryCif(Box::new(
            molframe_bcif::BcifBatchSource::new(
                source,
                options.clone(),
                DatasetId::new(0),
                ChunkId::new(0),
                LogicalRow::new(0),
                window_bytes,
                context,
            )?,
        ))),
        other => Err(StructureBatchError::Diagnostic(
            Diagnostic::new(Code::E1001)
                .with_message("format has no bounded structural batch reader")
                .with_context("format", other.name()),
        )),
    }
}

fn detect_format(
    source: &mut WindowedSourceFile,
    path: &Path,
    requested: Format,
    window_bytes: usize,
) -> Result<Format, StructureBatchError> {
    if requested != Format::Auto {
        return Ok(requested);
    }
    let window = source.window(0, window_bytes)?;
    for format in [Format::Mmcif, Format::BinaryCif, Format::Pdb] {
        if format.recognises(window.bytes()) {
            return Ok(format);
        }
    }
    let name = path.file_name().and_then(|name| name.to_str());
    name.and_then(Format::from_name).ok_or_else(|| {
        StructureBatchError::Diagnostic(
            Diagnostic::new(Code::E1001).with_context("name", path.display().to_string()),
        )
    })
}

fn source_window_bytes(context: &ExecutionContext) -> Result<usize, StructureBatchError> {
    let available = context
        .memory_budget()
        .bytes()
        .saturating_sub(context.reserved_bytes());
    let window = DEFAULT_SOURCE_WINDOW_BYTES.min(available / 4);
    if window < MINIMUM_SOURCE_WINDOW_BYTES {
        return Err(StructureBatchError::RecordExceedsBudget {
            required: MINIMUM_SOURCE_WINDOW_BYTES.saturating_mul(4),
            available,
        });
    }
    Ok(window)
}

#[cfg(test)]
#[path = "structure_batches_tests.rs"]
mod tests;
