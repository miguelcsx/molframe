//! Database sequence references and canonical alignment intervals.

use crate::{Category, DataBlock};
use pdbiox_core::structure::{ReferenceAlignment, ReferenceSequence, SequenceReferences};
use pdbiox_core::{Code, Diagnostic, Diagnostics};

pub(super) fn read(block: &DataBlock, findings: &mut Diagnostics) -> SequenceReferences {
    SequenceReferences {
        sequences: sequences(block, findings),
        alignments: alignments(block, findings),
    }
}

fn sequences(block: &DataBlock, findings: &mut Diagnostics) -> Vec<ReferenceSequence> {
    let Some(category) = block.category("struct_ref") else {
        return Vec::new();
    };
    let mut output = Vec::new();
    for row in 0..category.row_count() {
        let (Some(id), Some(entity_id)) = (
            category.identifier("id", row),
            category.identifier("entity_id", row),
        ) else {
            findings.push(row_error(category, row));
            continue;
        };
        output.push(ReferenceSequence {
            id: id.into_owned().into_boxed_str(),
            entity_id: entity_id.into_owned().into_boxed_str(),
            database_name: boxed(category.identifier("db_name", row)),
            database_code: boxed(category.identifier("db_code", row)),
            accession: boxed(category.identifier("pdbx_db_accession", row)),
            one_letter_code: category
                .text("pdbx_seq_one_letter_code", row)
                .map(|sequence| {
                    sequence
                        .chars()
                        .filter(|character| !character.is_whitespace())
                        .collect::<String>()
                        .into_boxed_str()
                }),
        });
    }
    output
}

fn alignments(block: &DataBlock, findings: &mut Diagnostics) -> Vec<ReferenceAlignment> {
    let Some(category) = block.category("struct_ref_seq") else {
        return Vec::new();
    };
    let mut output = Vec::new();
    for row in 0..category.row_count() {
        let values = (
            category.identifier("align_id", row),
            category.identifier("ref_id", row),
            integer(category, "seq_align_beg", row),
            integer(category, "seq_align_end", row),
            integer(category, "db_align_beg", row),
            integer(category, "db_align_end", row),
        );
        let (Some(id), Some(reference_id), Some(begin), Some(end), Some(db_begin), Some(db_end)) =
            values
        else {
            findings.push(row_error(category, row));
            continue;
        };
        let chain_ids: Box<[Box<str>]> = match category.identifier("pdbx_strand_id", row) {
            Some(value) => value
                .split(',')
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(Box::<str>::from)
                .collect(),
            None => Box::new([]),
        };
        output.push(ReferenceAlignment {
            id: id.into_owned().into_boxed_str(),
            reference_id: reference_id.into_owned().into_boxed_str(),
            chain_ids,
            canonical: [begin, end],
            reference: [db_begin, db_end],
        });
    }
    output
}

fn integer(category: &Category, item: &str, row: usize) -> Option<i32> {
    category
        .value(item, row)?
        .as_integer()
        .and_then(|value| i32::try_from(value).ok())
}

fn boxed(value: Option<std::borrow::Cow<'_, str>>) -> Option<Box<str>> {
    value.map(|value| value.into_owned().into_boxed_str())
}

fn row_error(category: &Category, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2001)
        .in_category(category.name())
        .at_row(u32::try_from(row).map_or(u32::MAX, |value| value))
}

#[cfg(test)]
#[path = "references_tests.rs"]
mod tests;
