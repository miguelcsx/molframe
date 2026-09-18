//! The CIF family's buffered readers.
//!
//! `mmCIF` and `BinaryCIF` decode different bytes and take the same three modes:
//! atomic coordinates only, the document projected through modelcif, or the
//! document materialized as it stands. The body is written once here and
//! invoked with the crate that owns the decoding, rather than copied per format.
//! `Document` is the same type in both, so nothing else differs.

use molframe_core::diagnostic::{Diagnostic, Findings};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

use super::extensions;

/// Defines the buffered reader of one crate in the CIF family.
macro_rules! cif_family_reader {
    ($reader:ident, $family:ident) => {
        pub(super) fn $reader(
            input: &InputBuffer,
            options: &ReadOptions,
        ) -> Result<(Structure, Vec<Diagnostic>), Findings> {
            if options.only_atomic_coords {
                return $family::read(input, options).map_err(Findings::from);
            }
            #[cfg(feature = "modelcif")]
            {
                let projection = match molframe_modelcif::ModelCifProjection::new(
                    molframe_modelcif::ModelCifOptions::new(),
                ) {
                    Ok(projection) => projection,
                    Err(error) => return Err(extensions::model_error(&error).into()),
                };
                let (document, structure, projected, findings) = $family::read_with_projection(
                    input,
                    options,
                    extensions::keep_non_model_extension_category,
                    projection,
                )?;
                extensions::attach_projected_metadata(
                    &document, structure, findings, projected, options,
                )
                .map_err(Findings::from)
            }
            #[cfg(not(feature = "modelcif"))]
            {
                let (document, structure, findings) = $family::read_with_metadata(
                    input,
                    options,
                    extensions::keep_non_model_extension_category,
                )?;
                extensions::attach_materialized_cif_metadata(
                    &document, structure, findings, options,
                )
                .map_err(Findings::from)
            }
        }
    };
}

#[cfg(feature = "bcif")]
cif_family_reader!(read_bcif_buffer, molframe_bcif);
#[cfg(feature = "mmcif")]
cif_family_reader!(read_mmcif_buffer, molframe_cif);
