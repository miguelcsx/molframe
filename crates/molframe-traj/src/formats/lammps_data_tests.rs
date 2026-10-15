use super::{LammpsAtomStyle, LammpsDataError, parse_lammps_data};

const DATA: &str = "water topology\n\n3 atoms\n2 bonds\n1 angles\n\n0 10 xlo xhi\n0 11 ylo yhi\n0 12 zlo zhi\n1.0 0.5 -0.25 xy xz yz\n\nMasses\n\n1 15.9994 # oxygen\n2 1.008\n\nAtoms # full\n\n30 7 2 0.417 1 2 3 0 0 0\n10 7 1 -0.834 0 2 3 1 -1 0\n20 7 2 0.417 2 2 3 0 0 0\n\nBonds\n\n2 1 10 30\n1 1 10 20\n\nAngles\n\n1 1 20 10 30\n\nPair Coeffs\n\n1 0.1 3.0\n";

#[test]
fn full_style_preserves_sparse_ids_cell_force_field_data_and_unknown_sections() {
    let topology = parse_lammps_data(DATA).expect("valid full-style data should parse");
    assert_eq!(topology.atom_style, LammpsAtomStyle::Full);
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.id)
            .collect::<Vec<_>>(),
        [10, 20, 30]
    );
    assert_eq!(topology.atoms[0].image, Some([1, -1, 0]));
    assert_eq!(topology.masses.get(&1), Some(&15.9994));
    assert_eq!(topology.bonds[0].atoms, [10, 20]);
    let tilt = topology.cell.expect("cell should exist").tilt;
    for (observed, expected) in tilt.into_iter().zip([1.0, 0.5, -0.25]) {
        assert!((observed - expected).abs() < f64::EPSILON);
    }
    assert_eq!(
        topology.other_sections["Pair Coeffs"][0].as_ref(),
        "1 0.1 3.0"
    );
}

#[test]
fn absent_atom_style_is_refused_because_six_column_records_are_ambiguous() {
    let source = "ambiguous\n1 atoms\n\nAtoms\n\n1 1 0 0 0\n";
    assert_eq!(
        parse_lammps_data(source),
        Err(LammpsDataError::MissingAtomStyle)
    );
}

#[test]
fn connectivity_to_an_absent_atom_is_rejected() {
    let source = "bad bond\n1 atoms\n1 bonds\n\nAtoms # atomic\n\n1 1 0 0 0\n\nBonds\n\n1 1 1 2\n";
    assert_eq!(parse_lammps_data(source), Err(LammpsDataError::UnknownAtom));
}
