//! Canonical database sequence-reference output.

use super::value::quote_text;
use pdbiox_core::structure::{SEQUENCE_REFERENCES_EXTENSION, SequenceReferences, Structure};
use std::fmt::Write as _;

pub(super) fn write(out: &mut String, structure: &Structure) {
    let Some(references) = structure
        .extensions()
        .get::<SequenceReferences>(SEQUENCE_REFERENCES_EXTENSION)
    else {
        return;
    };
    write_sequences(out, references);
    write_alignments(out, references);
}

fn write_sequences(out: &mut String, references: &SequenceReferences) {
    if references.sequences.is_empty() {
        return;
    }
    out.push_str(
        "loop_\n\
_struct_ref.id\n\
_struct_ref.entity_id\n\
_struct_ref.db_name\n\
_struct_ref.db_code\n\
_struct_ref.pdbx_db_accession\n\
_struct_ref.pdbx_seq_one_letter_code\n",
    );
    for sequence in &references.sequences {
        let _ = writeln!(
            out,
            "{} {} {} {} {} {}",
            quote_text(&sequence.id),
            quote_text(&sequence.entity_id),
            optional(sequence.database_name.as_deref()),
            optional(sequence.database_code.as_deref()),
            optional(sequence.accession.as_deref()),
            optional(sequence.one_letter_code.as_deref()),
        );
    }
    out.push_str("#\n");
}

fn write_alignments(out: &mut String, references: &SequenceReferences) {
    if references.alignments.is_empty() {
        return;
    }
    out.push_str(
        "loop_\n\
_struct_ref_seq.align_id\n\
_struct_ref_seq.ref_id\n\
_struct_ref_seq.pdbx_strand_id\n\
_struct_ref_seq.seq_align_beg\n\
_struct_ref_seq.seq_align_end\n\
_struct_ref_seq.db_align_beg\n\
_struct_ref_seq.db_align_end\n",
    );
    for alignment in &references.alignments {
        let chain_ids = alignment.chain_ids.join(",");
        let _ = writeln!(
            out,
            "{} {} {} {} {} {} {}",
            quote_text(&alignment.id),
            quote_text(&alignment.reference_id),
            optional((!chain_ids.is_empty()).then_some(chain_ids.as_str())),
            alignment.canonical[0],
            alignment.canonical[1],
            alignment.reference[0],
            alignment.reference[1],
        );
    }
    out.push_str("#\n");
}

fn optional(value: Option<&str>) -> String {
    match value {
        Some(value) => quote_text(value),
        None => "?".to_owned(),
    }
}
