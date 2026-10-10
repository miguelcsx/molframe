use super::*;
use molframe_core::io::{InputBuffer, ReadOptions};
use std::fmt::Write;

fn structure(second_contact: bool, coplanar: bool) -> Structure {
    let mut source = String::from(
        "data_pairs\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let atoms = [
        ("C1'", "DA", "A", [-5.0, 0.0, 0.0]),
        ("C2", "DA", "A", [0.0, 0.0, 0.0]),
        ("C4", "DA", "A", [0.0, 1.0, 0.0]),
        ("C6", "DA", "A", [1.0, 0.0, 0.0]),
        ("N6", "DA", "A", [0.0, 0.0, 0.0]),
        ("N1", "DA", "A", [0.0, 2.0, 0.0]),
        ("C1'", "DT", "B", [5.0, 0.0, 0.0]),
        ("C2", "DT", "B", [2.8, 0.0, 0.0]),
        (
            "C4",
            "DT",
            "B",
            if coplanar {
                [2.8, 1.0, 0.0]
            } else {
                [2.8, 0.0, 1.0]
            },
        ),
        ("C6", "DT", "B", [3.8, 0.0, 0.0]),
        ("O4", "DT", "B", [2.8, 0.0, 0.0]),
        (
            "N3",
            "DT",
            "B",
            [if second_contact { 2.8 } else { 20.0 }, 2.0, 0.0],
        ),
    ];
    for (row, (name, residue, chain, [x, y, z])) in atoms.into_iter().enumerate() {
        writeln!(
            source,
            "ATOM {} C {name} {residue} {chain} 1 {x} {y} {z}",
            row + 1
        )
        .expect("source writes");
    }
    molframe_cif::read(
        &InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("fixture reads")
    .0
}

#[test]
fn complementary_coplanar_bases_need_two_canonical_edge_contacts() {
    let pairs =
        standard_base_pair_geometry(&structure(true, true), WatsonCrickGeometry::default(), 3.6)
            .expect("geometry resolves");
    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0].contact_count, 2);
    assert!((pairs[0].closest_distance - 2.8).abs() < 1e-5);
    for structure in [structure(false, true), structure(true, false)] {
        assert!(
            standard_base_pair_geometry(&structure, WatsonCrickGeometry::default(), 3.6)
                .expect("geometry resolves")
                .is_empty()
        );
    }
}

#[test]
fn invalid_pairing_limits_return_a_spatial_error() {
    assert!(
        standard_base_pair_geometry(
            &structure(true, true),
            WatsonCrickGeometry::default(),
            f32::NAN
        )
        .is_err()
    );
}
