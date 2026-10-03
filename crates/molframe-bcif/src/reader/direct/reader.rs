//! One-pass direct `BinaryCIF` read entry points.

use super::container;
use super::ensemble::classify;
use super::projection::{self, decode, decode_with_projection};
use molframe_cif::{
    CifEventSink, Document, lower_atom_site_with, lower_ragged_atom_site_with,
    lower_single_atom_site_with,
};
use molframe_core::diagnostic::Diagnostic;
use molframe_core::io::{InputBuffer, ReadOptions, ReadResult};
use molframe_core::structure::Structure;

type ProjectionResult<T> = Result<(Document, Structure, T, Vec<Diagnostic>), Vec<Diagnostic>>;

pub(in crate::reader) fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    read_with_metadata(input, options, |_| false)
        .map(|(_, structure, findings)| (structure, findings))
}

pub(in crate::reader) fn read_with_metadata(
    input: &InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
) -> Result<(Document, Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
    decode_columns_with(input, options, keep_category)?.lower(options)
}

/// Container parse and column decode, held until lowering.
///
/// Splits [`read`] in two so the benches can time the stages separately.
#[doc(hidden)]
pub struct DecodedColumns<'a> {
    projection: projection::Projection<'a>,
}

impl std::fmt::Debug for DecodedColumns<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("DecodedColumns")
            .finish_non_exhaustive()
    }
}

/// Decodes the container and atom columns without lowering them.
///
/// # Errors
///
/// Returns ordered container or codec diagnostics.
#[doc(hidden)]
pub fn decode_columns<'a>(
    input: &'a InputBuffer,
    options: &ReadOptions,
) -> Result<DecodedColumns<'a>, Vec<Diagnostic>> {
    decode_columns_with(input, options, |_| false)
}

fn decode_columns_with<'a>(
    input: &'a InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
) -> Result<DecodedColumns<'a>, Vec<Diagnostic>> {
    let binary =
        container::parse(input.as_bytes(), options.limits).map_err(|finding| vec![finding])?;
    let projection = decode(binary, options, keep_category).map_err(|finding| vec![finding])?;
    Ok(DecodedColumns { projection })
}

impl DecodedColumns<'_> {
    /// Lowers the decoded columns into a structure.
    ///
    /// # Errors
    ///
    /// Returns ordered interpretation diagnostics.
    #[doc(hidden)]
    pub fn lower(
        self,
        options: &ReadOptions,
    ) -> Result<(Document, Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
        finish_projection(self.projection, options)
    }
}

pub(in crate::reader) fn read_with_projection<S>(
    input: &InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
    sink: S,
) -> ProjectionResult<S::Output>
where
    S: CifEventSink,
{
    let binary =
        container::parse(input.as_bytes(), options.limits).map_err(|finding| vec![finding])?;
    let (projection, projected) = decode_with_projection(binary, options, keep_category, sink)
        .map_err(|finding| vec![finding])?;
    let (metadata, structure, findings) = finish_projection(projection, options)?;
    Ok((metadata, structure, projected, findings))
}

fn finish_projection(
    projection: projection::Projection<'_>,
    options: &ReadOptions,
) -> Result<(Document, Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
    let projection::Projection {
        metadata,
        mut atoms,
        extra_atom_columns,
    } = projection;
    let layout = if options.only_first_model {
        None
    } else {
        Some(classify(&atoms, extra_atom_columns).map_err(|finding| vec![finding])?)
    };
    atoms.compact_floats_for_lowering();
    let lowered = match layout {
        Some(layout) if layout.ragged => lower_ragged_atom_site_with(
            &metadata,
            options,
            Vec::new(),
            layout.model_count,
            |sink| {
                atoms.feed(sink);
                Ok(())
            },
        ),
        None => {
            lower_single_atom_site_with(&metadata, options, Vec::new(), atoms.row_count(), |sink| {
                atoms.feed(sink);
                Ok(())
            })
        }
        Some(layout) if layout.model_count <= 1 => {
            lower_single_atom_site_with(&metadata, options, Vec::new(), atoms.row_count(), |sink| {
                atoms.feed(sink);
                Ok(())
            })
        }
        _ => lower_atom_site_with(&metadata, options, Vec::new(), |sink| {
            atoms.feed(sink);
            Ok(())
        }),
    };
    let (structure, findings) = lowered?;
    Ok((metadata, structure, findings))
}
