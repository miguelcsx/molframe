use super::{GromacsItpError, parse_gromacs_itp};

const ITP: &str = "#define FLEXIBLE\n[ moleculetype ]\nSOL 2\n\n[ atoms ]\n3 HW 1 SOL HW2 1 0.417 1.008\n1 OW 1 SOL OW 1 -0.834 15.9994\n2 HW 1 SOL HW1 1 0.417 1.008\n\n[ bonds ]\n1 2 1 0.09572 345000\n1 3 1 0.09572 345000\n\n[ angles ]\n2 1 3 1 1.824218 383\n\n[ exclusions ]\n1 2 3\n\n[ position_restraints ]\n1 1 1000 1000 1000\n";

#[test]
fn topology_preserves_directives_parameters_and_unknown_sections() {
    let topology = parse_gromacs_itp(ITP).expect("valid include topology should parse");
    assert_eq!(topology.directives[0].as_ref(), "#define FLEXIBLE");
    assert_eq!(
        topology
            .molecule_type
            .expect("molecule should exist")
            .name
            .as_ref(),
        "SOL"
    );
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.number)
            .collect::<Vec<_>>(),
        [1, 2, 3]
    );
    assert_eq!(topology.bonds[0].parameters[0].as_ref(), "0.09572");
    assert_eq!(topology.exclusions[0], [1, 2, 3]);
    assert_eq!(
        topology.other_sections["position_restraints"][0].as_ref(),
        "1 1 1000 1000 1000"
    );
}

#[test]
fn interaction_to_an_absent_atom_is_rejected() {
    let source = "[ atoms ]\n1 C 1 MOL C 1 0.0 12.0\n[ bonds ]\n1 2 1\n";
    assert_eq!(parse_gromacs_itp(source), Err(GromacsItpError::UnknownAtom));
}

#[test]
fn duplicate_atom_numbers_are_rejected() {
    let source = "[ atoms ]\n1 C 1 MOL C1 1 0.0\n1 C 1 MOL C2 1 0.0\n";
    assert_eq!(
        parse_gromacs_itp(source),
        Err(GromacsItpError::DuplicateAtom)
    );
}
