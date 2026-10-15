use super::{AmberTopologyError, parse_amber_topology};

const TOP: &str = "%VERSION  VERSION_STAMP = V0001.000\n%FLAG POINTERS\n%FORMAT(10I8)\n       3       2\n%FLAG ATOM_NAME\n%FORMAT(20a4)\nO   H1  H2  \n%FLAG CHARGE\n%FORMAT(5E16.8)\n -0.15191396E+02  0.75956980E+01  0.75956980E+01\n%FLAG MASS\n%FORMAT(5E16.8)\n  0.15999400E+02  0.10080000E+01  0.10080000E+01\n%FLAG ATOM_TYPE_INDEX\n%FORMAT(10I8)\n       1       2       2\n%FLAG AMBER_ATOM_TYPE\n%FORMAT(20a4)\nOW  HW  HW  \n%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n       8       1       1\n%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\nWAT \n%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n       1\n%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n       0       3       1       0       6       1\n%FLAG BONDS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n%FLAG ANGLES_INC_HYDROGEN\n%FORMAT(10I8)\n       3       0       6       1\n%FLAG ANGLES_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n%FLAG DIHEDRALS_INC_HYDROGEN\n%FORMAT(10I8)\n       0       3      -6      -3       1\n%FLAG DIHEDRALS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n%FLAG BOND_FORCE_CONSTANT\n%FORMAT(5E16.8)\n  0.55300000E+03\n";

#[test]
fn fixed_width_arrays_join_into_atoms_residues_and_signed_connectivity() {
    let topology = parse_amber_topology(TOP).expect("valid parm7 should parse");
    assert_eq!(topology.atoms.len(), 3);
    assert!((topology.atoms[0].charge + 0.833_670_61).abs() < 1e-7);
    assert_eq!(topology.atoms[0].atomic_number, Some(8));
    assert_eq!(topology.residues[0].atom_start, 0);
    assert_eq!(topology.residues[0].atom_end, 3);
    assert_eq!(topology.bonds[1].atoms, [0, 2]);
    assert_eq!(topology.angles[0].atoms, [1, 0, 2]);
    assert!(topology.dihedrals[0].improper);
    assert!(topology.dihedrals[0].ignore_end_group);
    assert_eq!(
        topology.sections["BOND_FORCE_CONSTANT"].values[0].as_ref(),
        "0.55300000E+03"
    );
}

#[test]
fn pointer_count_must_match_parallel_atom_arrays() {
    let broken = TOP.replacen("       3       2", "       4       2", 1);
    assert_eq!(
        parse_amber_topology(&broken),
        Err(AmberTopologyError::LengthMismatch)
    );
}

#[test]
fn coordinate_offsets_must_be_multiples_of_three() {
    let broken = TOP.replacen("       0       3       1", "       0       4       1", 1);
    assert_eq!(
        parse_amber_topology(&broken),
        Err(AmberTopologyError::InvalidAtomReference)
    );
}

#[test]
fn external_parmed_fixture_when_configured() {
    let Ok(path) = std::env::var("MOLFRAME_AMBER_TOP_FIXTURE") else {
        return;
    };
    let source = std::fs::read_to_string(path).expect("external parm7 fixture should be readable");
    let topology = parse_amber_topology(&source).expect("external parm7 fixture should parse");
    assert!(!topology.atoms.is_empty());
    assert!(!topology.residues.is_empty());
    assert!(topology.atoms.iter().all(|atom| {
        usize::try_from(atom.residue).is_ok_and(|residue| residue < topology.residues.len())
    }));
}
