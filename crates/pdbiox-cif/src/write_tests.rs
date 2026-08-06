use super::*;
use crate::document::{Category, DataBlock};
use crate::lexer::Quoting;
use crate::{Document, parse, read};
use pdbiox_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "\
data_TEST
_entry.id   TEST
#
_my_lab.note   'kept verbatim'
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N  . GLY A 1 ? 27.340 24.430 2.614 1.00 10.00 1 A 1
ATOM 2 C CA . GLY A 1 ? 26.266 25.413 2.842 1.00 11.00 1 A 1
ATOM 3 N N  . ASN A 2 ? 26.335 27.770 3.258 1.00 14.00 2 A 1
#
";

fn structure(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

fn document(text: &str) -> Document {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    }
}

#[test]
fn a_canonical_write_reads_back_to_the_same_structure() {
    let original = structure(SOURCE);
    let round_tripped = structure(&write_canonical(&original));

    assert_eq!(round_tripped.atom_count(), original.atom_count());
    assert_eq!(round_tripped.residue_count(), original.residue_count());
    assert_eq!(round_tripped.chain_count(), original.chain_count());
}

#[test]
fn positions_survive_a_canonical_round_trip_to_the_precision_written() {
    let original = structure(SOURCE);
    let round_tripped = structure(&write_canonical(&original));
    for (before, after) in original.positions().iter().zip(round_tripped.positions()) {
        for axis in 0..3 {
            assert!(
                (before[axis] - after[axis]).abs() < 1e-3,
                "{before:?} became {after:?}"
            );
        }
    }
}

#[test]
fn a_preserving_write_keeps_a_category_the_library_does_not_interpret() {
    let written = write_preserving(&document(SOURCE));
    assert!(written.contains("_my_lab.note"), "{written}");
    assert!(written.contains("kept verbatim"));
}

#[test]
fn a_preserving_write_reparses_to_the_same_document() {
    let original = document(SOURCE);
    let round_tripped = document(&write_preserving(&original));

    let items = |doc: &Document| -> Vec<String> {
        doc.first_block()
            .into_iter()
            .flat_map(DataBlock::categories)
            .flat_map(|category: &Category| {
                let name = category.name().to_owned();
                category
                    .items()
                    .map(move |item| format!("{name}.{item}"))
                    .collect::<Vec<_>>()
            })
            .collect()
    };
    assert_eq!(items(&round_tripped), items(&original));
}

#[test]
fn a_value_that_would_be_reread_as_something_else_is_quoted() {
    let mut block = DataBlock::new("x");
    let category = block.category_mut("a", pdbiox_core::span::ByteSpan::default());
    category
        .column_mut("spaced")
        .push(CifValue::Text("two words".into()), Quoting::Bare);
    let mut document = Document::new();
    document.push(block);

    let written = write_preserving(&document);
    assert!(written.contains("'two words'"), "{written}");
    let reparsed = document_value(&written);
    assert_eq!(reparsed.as_deref(), Some("two words"));
}

fn document_value(text: &str) -> Option<String> {
    document(text)
        .first_block()?
        .category("a")?
        .text("spaced", 0)
        .map(str::to_owned)
}

#[test]
fn a_sentinel_survives_a_preserving_round_trip_as_a_sentinel() {
    let written = write_preserving(&document("data_x\n_a.gone .\n_a.unknown ?\n"));
    let back = document(&written);
    let category = back.first_block().and_then(|block| block.category("a"));
    assert_eq!(
        category.and_then(|c| c.value("gone", 0)),
        Some(&CifValue::Inapplicable)
    );
    assert_eq!(
        category.and_then(|c| c.value("unknown", 0)),
        Some(&CifValue::Unknown)
    );
}
