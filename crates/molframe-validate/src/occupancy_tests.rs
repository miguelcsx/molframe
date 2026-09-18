use super::*;
use molframe_core::io::{InputBuffer, ReadOptions};

fn structure() -> Structure {
    let cif = b"data_o
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
ATOM 1 C CA A ALA A 1 0 0 0 0.6
ATOM 2 C CA B ALA A 1 0 0 0 0.4
ATOM 3 C CB A ALA A 1 1 0 0 0.7
ATOM 4 C CB B ALA A 1 1 0 0 0.2
";
    let input = InputBuffer::from_bytes(cif.to_vec());
    molframe_cif::read(&input, &ReadOptions::new()).map_or_else(
        |findings| panic!("fixture failed: {findings:?}"),
        |result| result.0,
    )
}

#[test]
fn alternate_atom_sums_are_grouped_without_component_tables() {
    let report = altloc_occupancy_sums(
        &structure(),
        Namespace::Label,
        AltlocOccupancyOptions {
            expected_sum: 1.0,
            tolerance: 1e-6,
        },
    )
    .expect("valid identifiers");
    assert_eq!(report.intended, 4);
    assert_eq!(report.assessed, 4);
    assert_eq!(report.records.len(), 2);
    assert_eq!(report.records[0].issue, None);
    assert_eq!(
        report.records[1].issue,
        Some(AltlocOccupancyIssue::SumMismatch)
    );
}
