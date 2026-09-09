//! One-pass direct-read entry points and layout diagnostics.

use super::extra::ProjectedSink;
use super::projection;
use super::stream::DirectSink;
use crate::document::Document;
use crate::lower::finish_streamed;
use crate::parser::{CifEventSink, parse_into};
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::io::{InputBuffer, ReadOptions, ReadResult};
use pdbiox_core::structure::Structure;

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
    let (direct, findings) =
        parse_into(input, DirectSink::new(options, keep_category, input.len()))?;
    finish_direct(direct, findings, options)
}

pub(in crate::reader) fn read_with_projection<S>(
    input: &InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
    projection: S,
) -> ProjectionResult<S::Output>
where
    S: CifEventSink,
{
    let sink = ProjectedSink::new(options, keep_category, input.len(), projection);
    let ((direct, projected), findings) = parse_into(input, sink)?;
    let (document, structure, findings) = finish_direct(direct, findings, options)?;
    Ok((document, structure, projected, findings))
}

fn finish_direct(
    direct: super::stream::DirectOutput,
    findings: Vec<Diagnostic>,
    options: &ReadOptions,
) -> Result<(Document, Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
    let projection = direct.projection;
    if !projection.layout.has_atom_site {
        return Err(vec![
            Diagnostic::new(Code::E2001)
                .with_message("the block contains no coordinates")
                .in_category("atom_site"),
        ]);
    }
    if !projection.layout.supports_direct() {
        return Err(unsupported_layout(&projection.layout, findings));
    }

    if !direct.setup_errors.is_empty() {
        return Err(direct.setup_errors);
    }
    let Some(models) = direct.models else {
        return Err(vec![
            Diagnostic::new(Code::E2001)
                .with_message("the block contains no coordinates")
                .in_category("atom_site"),
        ]);
    };
    let (structure, findings) = finish_streamed(&projection.metadata, options, findings, models)?;
    Ok((projection.metadata, structure, findings))
}

fn unsupported_layout(
    layout: &projection::AtomLayout,
    findings: Vec<Diagnostic>,
) -> Vec<Diagnostic> {
    if findings.iter().any(|finding| finding.code() == Code::E1103) {
        return findings;
    }
    let message = if layout.model_count() == 0 {
        "atom_site does not contain complete coordinate rows"
    } else {
        "atom_site is split across scalar, repeated, or duplicate loop items"
    };
    let mut diagnostics = Diagnostics::with_capacity(findings.len() + 1);
    diagnostics.extend(findings);
    diagnostics.push(
        Diagnostic::new(Code::E1104)
            .with_message(message)
            .in_category("atom_site"),
    );
    diagnostics.finish()
}

pub(in crate::reader) fn keep_lowering_category(category: &str) -> bool {
    matches!(
        category,
        "entry"
            | "struct"
            | "refine"
            | "exptl"
            | "cell"
            | "entity"
            | "entity_poly_seq"
            | "struct_asym"
            | "struct_conn"
            | "struct_ref"
            | "struct_ref_seq"
    )
}
