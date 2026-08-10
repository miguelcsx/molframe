use super::{AimsError, parse_aims_geometry, write_aims_geometry};

#[test]
fn fractional_atoms_use_all_three_skewed_lattice_vectors() {
    let source = "lattice_vector 2 0 0\n\
                  lattice_vector 1 3 0\n\
                  lattice_vector 0.5 0.25 4\n\
                  atom_frac 0.5 0.5 0.5 C\n";
    let geometry = parse_aims_geometry(source)
        .unwrap_or_else(|error| panic!("FHI-aims geometry failed: {error}"));
    assert!(
        geometry.atoms[0]
            .position
            .iter()
            .zip([1.75, 1.625, 2.0])
            .all(|(observed, expected)| (observed - expected).abs() < f32::EPSILON)
    );
    let cell = geometry.cell().unwrap_or_else(|| panic!("cell absent"));
    assert!((cell.lengths[1] - 10.0_f64.sqrt()).abs() < 1.0e-12);
}

#[test]
fn writer_round_trips_species_coordinates_and_cell() {
    let source = "# molecule\nlattice_vector 10 0 0\n\
                  lattice_vector 0 11 0\nlattice_vector 0 0 12\n\
                  atom -1.25 2.5 3.75 Cl\natom 4 5 6 H\n";
    let expected = parse_aims_geometry(source)
        .unwrap_or_else(|error| panic!("FHI-aims geometry failed: {error}"));
    let observed = parse_aims_geometry(&write_aims_geometry(&expected))
        .unwrap_or_else(|error| panic!("written FHI-aims geometry failed: {error}"));
    assert_eq!(observed, expected);
}

#[test]
fn fractional_atom_without_complete_cell_is_rejected() {
    let result = parse_aims_geometry("lattice_vector 1 0 0\natom_frac 0 0 0 C\n");
    assert_eq!(result, Err(AimsError::FractionalWithoutCell));
}
