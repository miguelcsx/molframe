//! First-block metadata projection and selected coordinate decoding.

use super::column::{AtomColumns, field_index};
use super::container::{DirectCategory, DirectColumn, DirectFile, decode_category};
use super::external::{feed_decoded, feed_encoded};
use pdbiox_cif::{CifEventSink, CifScalar, DataBlock, Document};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::span::ByteSpan;
use std::collections::HashSet;

pub(super) struct Projection<'a> {
    pub(super) metadata: Document,
    pub(super) atoms: AtomColumns,
    pub(super) extra_atom_columns: Vec<DirectColumn<'a>>,
}

pub(super) fn decode(
    binary: DirectFile<'_>,
    keep_category: fn(&str) -> bool,
) -> Result<Projection<'_>, Diagnostic> {
    let (projection, ()) = decode_with_projection(binary, keep_category, IgnoreSink)?;
    Ok(projection)
}

pub(super) fn decode_with_projection<S>(
    binary: DirectFile<'_>,
    keep_category: fn(&str) -> bool,
    mut sink: S,
) -> Result<(Projection<'_>, S::Output), Diagnostic>
where
    S: CifEventSink,
{
    let Some(encoded_block) = binary.data_blocks.into_iter().next() else {
        return Err(Diagnostic::new(Code::E1106));
    };
    sink.block(&encoded_block.header);
    let mut metadata_block = DataBlock::new(encoded_block.header);
    let mut atom_site = None;
    for encoded in encoded_block.categories {
        let name = encoded.name.trim_start_matches('_');
        let keep = keep_lowering_category(name) || keep_category(name);
        let project = sink.accepts_category(name);
        if name == "atom_site" {
            if atom_site.is_some() {
                return Err(Diagnostic::new(Code::E1401)
                    .with_message("a BinaryCIF block contains duplicate atom_site categories"));
            }
            atom_site = Some(encoded);
        } else if keep {
            let decoded = decode_category(encoded)?;
            if project {
                feed_decoded(&decoded, &mut sink);
            }
            let name = decoded.name().to_owned();
            *metadata_block.category_mut(&name, ByteSpan::default()) = decoded;
        } else if project {
            feed_encoded(encoded, &mut sink)?;
        }
    }
    let Some(atom_site) = atom_site else {
        return Err(Diagnostic::new(Code::E2001)
            .with_message("the block contains no coordinates")
            .in_category("atom_site"));
    };
    let (atoms, extra_atom_columns) = decode_atom_site(atom_site)?;
    let output = sink.finish();
    let mut metadata = Document::new();
    metadata.push(metadata_block);
    Ok((
        Projection {
            metadata,
            atoms,
            extra_atom_columns,
        },
        output,
    ))
}

struct IgnoreSink;

impl CifEventSink for IgnoreSink {
    type Output = ();

    fn block(&mut self, _name: &str) {}

    fn accepts_category(&self, _category: &str) -> bool {
        false
    }

    fn value(&mut self, _category: &str, _item: &str, _value: CifScalar<'_>, _span: ByteSpan) {}

    fn finish(self) {}
}

fn decode_atom_site(
    encoded: DirectCategory<'_>,
) -> Result<(AtomColumns, Vec<DirectColumn<'_>>), Diagnostic> {
    let mut atoms = AtomColumns::new(encoded.row_count);
    let mut extras = Vec::new();
    let mut names = HashSet::with_capacity(encoded.columns.len());
    for column in encoded.columns {
        if !names.insert(column.name.clone()) {
            return Err(Diagnostic::new(Code::E1401)
                .with_message("a BinaryCIF category contains a duplicate item")
                .with_context("item", column.name));
        }
        if field_index(&column.name).is_some() {
            let name = column.name.clone();
            atoms.insert(&name, column)?;
        } else {
            extras.push(column);
        }
    }
    Ok((atoms, extras))
}

fn keep_lowering_category(category: &str) -> bool {
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
