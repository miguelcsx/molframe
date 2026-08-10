use super::{QualityIssue, quality_flags};
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;

// Atom 0 is clean; atom 1 has zero occupancy, atom 2 occupancy above one, atom 3
// a negative B-factor.
const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n\
ATOM 1 C CA GLY A 1 0 0 0 1.0 20.0\n\
ATOM 2 C CB GLY A 2 1 0 0 0.0 20.0\n\
ATOM 3 C CG GLY A 3 2 0 0 1.5 20.0\n\
ATOM 4 C CD GLY A 4 3 0 0 1.0 -5.0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn each_out_of_range_value_is_flagged_against_its_atom() {
    let flags = quality_flags(&structure());
    let pairs: Vec<(u32, QualityIssue)> = flags
        .iter()
        .map(|flag| (flag.atom.get(), flag.issue))
        .collect();
    assert_eq!(
        pairs,
        vec![
            (1, QualityIssue::ZeroOccupancy),
            (2, QualityIssue::OccupancyOutOfRange),
            (3, QualityIssue::NegativeBFactor),
        ]
    );
}

#[test]
fn a_well_formed_atom_raises_no_flag() {
    let flags = quality_flags(&structure());
    assert!(flags.iter().all(|flag| flag.atom.get() != 0));
}
