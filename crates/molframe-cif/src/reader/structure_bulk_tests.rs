//! Large `atom_site` loops, which the reader takes in parallel pieces, must
//! read exactly as the serial parse does.

use super::*;

const HEADER: &str = "\
data_BULK
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
";

fn row(serial: usize, model: usize) -> String {
    let residue = serial / 4 + 1;
    let name = ["N", "CA", "C", "O"][serial % 4];
    let x = serial % 997;
    format!("ATOM {serial} C {name} GLY A {residue} {x}.5 {model}.25 -{x}.125 {model}\n")
}

/// Accepts every category, so the reader must deliver each `atom_site` value
/// to it and cannot take rows in bulk.
struct Everything;

impl crate::CifEventSink for Everything {
    type Output = ();

    fn block(&mut self, _name: &str) {}

    fn value(
        &mut self,
        _category: &str,
        _item: &str,
        _value: crate::CifScalar<'_>,
        _span: molframe_core::span::ByteSpan,
    ) {
    }

    fn finish(self) -> Self::Output {}
}

fn atoms(count: usize, models: usize) -> String {
    let mut text = String::from(HEADER);
    for model in 1..=models {
        for serial in 1..=count {
            text.push_str(&row(serial, model));
        }
    }
    text
}

fn assert_reads_like_the_serial_parse(text: &str) {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new();
    let bulk = read(&input, &options);
    let serial = read_with_projection(&input, &options, |_| false, Everything)
        .map(|(_, structure, (), findings)| (structure, findings));
    match (bulk, serial) {
        (Ok((bulk, bulk_findings)), Ok((serial, serial_findings))) => {
            let difference = structure_difference(
                &bulk,
                &serial,
                StructureDifferenceOptions {
                    coordinate_tolerance: 0.0,
                },
            )
            .expect("zero is a valid coordinate tolerance");
            assert!(difference.is_empty(), "difference: {difference:?}");
            assert_eq!(bulk.atom_count(), serial.atom_count());
            assert_eq!(bulk_findings, serial_findings);
        }
        (Err(bulk), Err(serial)) => assert_eq!(bulk, serial),
        (bulk, serial) => panic!(
            "one read failed: bulk ok {}, serial ok {}",
            bulk.is_ok(),
            serial.is_ok()
        ),
    }
}

#[test]
fn a_loop_cut_into_many_pieces_reads_in_order() {
    let text = atoms(9000, 2);
    assert!(text.len() > 1 << 19);
    assert_reads_like_the_serial_parse(&text);
}

#[test]
fn rows_spread_over_several_lines_read_whole() {
    let text = atoms(9000, 1).replace(" GLY A ", "\nGLY\nA ");
    assert_reads_like_the_serial_parse(&text);
}

#[test]
fn a_text_block_ends_the_parallel_read_where_it_starts() {
    let mut text = atoms(6000, 1);
    text = text.replacen("ATOM 3000 C N GLY", "ATOM 3000 C N\n;GLY\n;", 1);
    assert_reads_like_the_serial_parse(&text);
}

#[test]
fn a_following_loop_ends_the_parallel_read_at_its_header() {
    let mut text = atoms(9000, 1);
    text.push_str("loop_\n_other.a\n_other.b\n1 2\n3 4\n");
    assert_reads_like_the_serial_parse(&text);
}

#[test]
fn an_incomplete_last_row_is_reported_as_the_serial_parse_does() {
    let mut text = atoms(9000, 1);
    text.push_str("ATOM 9001 C N GLY\n");
    assert_reads_like_the_serial_parse(&text);
}
