use super::*;
use crate::{BinaryDocument, Decoded, decode, write_document};
use molframe_cif::{CifWriteOptions, parse, write_canonical_with_options};
use molframe_core::io::{InputBuffer, Limits, ReadOptions};
use molframe_core::structure::{StructureDifferenceOptions, structure_difference};
use std::fmt::Write as _;
use std::io::{self, Write};

const SIMPLE: &str = "\
data_DEMO
_entry.id DEMO
loop_
_entity.id
_entity.type
E1 polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N GLY A E1 1 1.125 2.25 3.5 1.0 12.5 7 GLY A N 1
ATOM 2 C CA GLY A E1 1 2.125 3.25 4.5 1.0 13.5 7 GLY A CA 1
";

#[test]
fn direct_structure_bytes_match_the_previous_projection_for_text_identifiers() {
    let structure = read_structure(SIMPLE);
    let options = CifWriteOptions::new();
    let actual = write_structure_with_options(&structure, &options).expect("direct write succeeds");
    let canonical = write_canonical_with_options(&structure, &options).expect("CIF writes");
    let input = InputBuffer::from_bytes(canonical.into_bytes());
    let (document, findings) = parse(&input).expect("canonical CIF parses");
    assert!(findings.is_empty(), "findings: {findings:?}");
    let expected = write_document(&document).expect("document writes");
    assert_eq!(actual, expected);
}

#[test]
fn numeric_and_text_identifiers_share_a_string_column_and_round_trip() {
    let source = SIMPLE.replace(
        "ATOM 2 C CA GLY A E1 1 2.125 3.25 4.5 1.0 13.5 7 GLY A CA 1",
        "ATOM 2 C CA GLY 1 E1 1 2.125 3.25 4.5 1.0 13.5 7 GLY 1 CA 1",
    );
    let structure = read_structure(&source);
    let first = write_structure(&structure).expect("mixed identifier column writes");
    let second = write_structure(&structure).expect("repeated write succeeds");
    assert_eq!(first, second);

    let input = InputBuffer::from_bytes(first);
    let (decoded, findings) = crate::read(&input, &ReadOptions::new()).expect("BCIF reads");
    assert!(findings.is_empty(), "findings: {findings:?}");
    let difference = structure_difference(
        &structure,
        &decoded,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("tolerance is valid");
    assert!(difference.is_empty(), "difference: {difference:?}");
}

#[derive(Default)]
struct CountingWriter {
    bytes: usize,
    writes: usize,
    largest_write: usize,
    checksum: u64,
}

impl Write for CountingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes += bytes.len();
        self.writes += 1;
        self.largest_write = self.largest_write.max(bytes.len());
        for byte in bytes {
            self.checksum = self.checksum.wrapping_mul(16_777_619) ^ u64::from(*byte);
        }
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn binary_writer_forwards_columns_without_a_complete_file_buffer() {
    let structure = read_structure(SIMPLE);
    let expected = write_structure(&structure).expect("in-memory write succeeds");
    let expected_checksum = expected.iter().fold(0_u64, |checksum, byte| {
        checksum.wrapping_mul(16_777_619) ^ u64::from(*byte)
    });
    let mut writer = CountingWriter::default();
    write_structure_to(&structure, &CifWriteOptions::new(), &mut writer)
        .expect("stream write succeeds");
    assert_eq!(writer.bytes, expected.len());
    assert_eq!(writer.checksum, expected_checksum);
    assert!(writer.writes > 1);
    assert!(writer.largest_write < expected.len());
}

#[test]
fn binary_workspace_limit_refuses_before_the_first_write() {
    let structure = read_structure(SIMPLE);
    let mut writer = CountingWriter::default();
    let result =
        write_structure_to_with_memory_limit(&structure, &CifWriteOptions::new(), 79, &mut writer);
    assert!(result.is_err());
    assert_eq!(writer.bytes, 0);
    assert_eq!(writer.writes, 0);
}

#[test]
fn a_large_repetitive_atom_column_writes_one_dictionary_value() {
    const ROWS: usize = 100_000;
    let mut source = String::from(
        "data_STRESS\n_entry.id STRESS\nloop_\n_entity.id\n_entity.type\nE1 polymer\n\
         loop_\n\
         _atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n\
         _atom_site.Cartn_x\n_atom_site.Cartn_y\n\
         _atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n",
    );
    for row in 1..=ROWS {
        let _ = writeln!(source, "ATOM {row} C CA GLY A E1 {row} 0 0 0 1");
    }
    let structure = read_structure(&source);
    let bytes = write_structure(&structure).expect("large structure writes");
    let binary = BinaryDocument::parse(&bytes, Limits::default()).expect("container parses");
    let atom_site = binary.file.data_blocks[0]
        .categories
        .iter()
        .find(|category| category.name == "_atom_site")
        .expect("atom-site category exists");
    let atom_name = atom_site
        .columns
        .iter()
        .find(|column| column.name == "label_atom_id")
        .expect("atom-name column exists");
    let Decoded::Strings(values) = decode(&atom_name.data).expect("column decodes") else {
        panic!("atom-name column is text")
    };
    assert_eq!(values.len(), ROWS);
    assert_eq!(values.dictionary().len(), 1);
    assert_eq!(values.dictionary()[0].as_ref(), "CA");
}

fn read_structure(source: &str) -> molframe_core::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, findings)) => {
            assert!(findings.is_empty(), "findings: {findings:?}");
            structure
        }
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}
