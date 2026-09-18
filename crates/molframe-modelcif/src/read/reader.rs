//! Public direct compact-read entry points.

use super::ModelCifProjection;
use crate::{ModelCif, ModelCifOptions, ModelCifReadError};
use molframe_cif::parse_events;
use molframe_core::{Diagnostic, io::InputBuffer};

/// Reads compact `ModelCIF` categories without constructing a CIF document.
///
/// Text is retained once per distinct dictionary value. Numeric values and
/// dictionary indices widen in place only when their observed range requires
/// it, so the retained representation remains bounded by the default policy.
///
/// # Errors
///
/// Returns syntax findings or a compact-storage policy failure.
pub fn read_compact(input: &InputBuffer) -> Result<(ModelCif, Vec<Diagnostic>), ModelCifReadError> {
    read_compact_with_options(input, ModelCifOptions::new())
}

/// Reads compact `ModelCIF` categories under an explicit memory policy.
///
/// # Errors
///
/// Returns syntax findings or a compact-storage policy failure.
pub fn read_compact_with_options(
    input: &InputBuffer,
    options: ModelCifOptions,
) -> Result<(ModelCif, Vec<Diagnostic>), ModelCifReadError> {
    let sink = ModelCifProjection::new(options).map_err(ModelCifReadError::Projection)?;
    let (projected, mut findings) = parse_events(input, sink).map_err(ModelCifReadError::Syntax)?;
    let (model, semantic) = projected.map_err(ModelCifReadError::Projection)?;
    findings.extend(semantic);
    Ok((model, findings))
}
