use super::*;
use crate::document::{Category, DataBlock};
use crate::lexer::Quoting;
use crate::{Document, parse, read};
use molframe_core::Structure;
use molframe_core::io::{InputBuffer, ReadOptions};
use proptest::prelude::*;

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

fn write_canonical(structure: &Structure) -> String {
    let options = crate::CifWriteOptions::new()
        .with_block_id("TEST")
        .with_generated_connection_ids()
        .with_connection_type_id("covale");
    match crate::write_canonical_with_options(structure, &options) {
        Ok(text) => text,
        Err(error) => panic!("canonical write failed: {error}"),
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
fn canonical_writing_does_not_fabricate_missing_author_identifiers() {
    let original = structure(SOURCE);
    let round_tripped = structure(&write_canonical(&original));
    let residue = round_tripped.residue(molframe_core::index::ResidueIndex::new(0));
    let atom = round_tripped.atom(molframe_core::index::AtomIndex::new(0));
    assert_eq!(
        residue.and_then(molframe_core::structure::ResidueRef::auth_name),
        None
    );
    assert_eq!(
        atom.and_then(molframe_core::structure::AtomRef::auth_name),
        None
    );
}

#[test]
fn canonical_writing_preserves_metadata_entities_and_both_name_spaces() {
    let source = "data_X\n\
        _entry.id X\n_struct.title 'A title'\n_exptl.method 'X-RAY DIFFRACTION'\n\
        _refine.ls_d_res_high 1.25\n\
        loop_\n_entity.id\n_entity.type\n_entity.pdbx_description\n7 polymer 'test entity'\n\
        loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n\
        7 1 GLY\n\
        loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
        _atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
        _atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
        _atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.auth_atom_id\n\
        _atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n\
        ATOM 1 C CA GLY LABEL 7 1 1 2 3 CAX GLYX AUTH 42\n";
    let original = structure(source);
    let round_tripped = structure(&write_canonical(&original));
    let entry = &round_tripped.data().entry;
    assert_eq!(entry.title.as_deref(), Some("A title"));
    assert_eq!(entry.method.as_deref(), Some("X-RAY DIFFRACTION"));
    assert!(
        entry
            .resolution
            .is_some_and(|value| (value - 1.25).abs() < f32::EPSILON)
    );
    let chain = round_tripped.chain(molframe_core::index::ChainIndex::new(0));
    assert_eq!(
        chain.and_then(molframe_core::structure::ChainRef::label),
        Some("LABEL")
    );
    assert_eq!(
        chain.and_then(molframe_core::structure::ChainRef::auth_label),
        Some("AUTH")
    );
    let residue = round_tripped.residue(molframe_core::index::ResidueIndex::new(0));
    assert_eq!(
        residue.and_then(molframe_core::structure::ResidueRef::auth_name),
        Some("GLYX")
    );
    let atom = round_tripped.atom(molframe_core::index::AtomIndex::new(0));
    assert_eq!(
        atom.and_then(molframe_core::structure::AtomRef::auth_name),
        Some("CAX")
    );
    assert_eq!(round_tripped.entity_count(), 1);
}

#[test]
fn a_preserving_write_keeps_a_category_the_library_does_not_interpret() {
    let written = write_preserving(&document(SOURCE));
    assert!(written.contains("_my_lab.note"), "{written}");
    assert!(written.contains("kept verbatim"));
}

#[test]
fn preserving_stream_matches_the_explicit_in_memory_wrapper() {
    let document = document(SOURCE);
    let expected = write_preserving(&document);
    let mut streamed = Vec::new();
    if let Err(error) = write_preserving_to(&document, &mut streamed) {
        panic!("preserving stream failed: {error}");
    }
    assert_eq!(streamed, expected.as_bytes());
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
    let category = block.category_mut("a", molframe_core::span::ByteSpan::default());
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

proptest! {
    #[test]
    fn generated_canonical_coordinates_survive_a_round_trip(
        positions in prop::collection::vec(
            (-9_000.0_f32..9_000.0, -9_000.0_f32..9_000.0, -9_000.0_f32..9_000.0),
            1..64,
        )
    ) {
        let source = generated_cif(&positions);
        let original = structure(&source);
        let round_tripped = structure(&write_canonical(&original));
        prop_assert_eq!(round_tripped.atom_count(), original.atom_count());
        prop_assert_eq!(round_tripped.residue_count(), original.residue_count());
        for (before, after) in original.positions().iter().zip(round_tripped.positions()) {
            for axis in 0..3 {
                // Canonical text records three decimal places, so half one unit
                // in the last place plus binary conversion error is the bound.
                prop_assert!((before[axis] - after[axis]).abs() <= 5.1e-4);
            }
        }
    }

    #[test]
    fn generated_document_values_survive_a_preserving_round_trip(
        value in "[A-Za-z][A-Za-z0-9_-]{0,40}"
    ) {
        let source = format!("data_generated\n_custom.value {value}\n");
        let original = document(&source);
        let round_tripped = document(&write_preserving(&original));
        let actual = round_tripped
            .first_block()
            .and_then(|block| block.category("custom"))
            .and_then(|category| category.text("value", 0));
        prop_assert_eq!(actual, Some(value.as_str()));
    }
}

fn generated_cif(positions: &[(f32, f32, f32)]) -> String {
    let mut source = String::from(
        "data_generated\nloop_\n\
         _atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n_atom_site.label_seq_id\n\
         _atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    for (index, (x, y, z)) in positions.iter().enumerate() {
        let _written = writeln!(
            source,
            "ATOM {} C C{} GLY A 1 {x:.3} {y:.3} {z:.3}",
            index + 1,
            index + 1,
        );
    }
    source.push_str("#\n");
    source
}
