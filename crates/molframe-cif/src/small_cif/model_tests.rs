use super::*;
use crate::parse;
use molframe_core::io::InputBuffer;

const CIF: &str = "data_small\n_cell_length_a 10\n_cell_length_b 11\n_cell_length_c 12\n_cell_angle_alpha 90\n_cell_angle_beta 91\n_cell_angle_gamma 92\n_symmetry_space_group_name_H-M 'P 1'\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\nC1 C 0.1 0.2 0.3 1\nO1 O 0.2 0.2 0.3 0.5\nloop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_distance\n_geom_bond_site_symmetry_2\nC1 O1 1.23 .\n";

fn parse_small(text: &str) -> Result<SmallCifStructure, SmallCifError> {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (document, findings) = parse(&input).expect("CIF syntax should parse");
    assert!(findings.is_empty());
    lower_small_cif_with_options(&document, SmallCifOptions::ddl1())
}

#[test]
fn ddl1_atoms_cell_fractional_coordinates_and_bonds_are_lowered() {
    let structure = parse_small(CIF).expect("core CIF should lower");
    assert_eq!(structure.atoms.len(), 2);
    assert_eq!(structure.atoms[1].fractional, Some([0.2, 0.2, 0.3]));
    assert_eq!(structure.atoms[1].occupancy, Some(0.5));
    assert_eq!(structure.bonds[0].first, 0);
    assert_eq!(structure.bonds[0].second, 1);
    assert!(
        structure
            .cell
            .expect("cell should exist")
            .lengths
            .into_iter()
            .zip([10.0, 11.0, 12.0])
            .all(|(actual, expected)| (actual - expected).abs() < f64::EPSILON)
    );
}

#[test]
fn ddl1_names_require_the_explicit_dialect() {
    let input = InputBuffer::from_bytes(CIF.as_bytes().to_vec());
    let (document, _) = parse(&input).expect("CIF syntax should parse");
    assert_eq!(lower_small_cif(&document), Err(SmallCifError::MissingAtoms));
}

#[test]
fn missing_atom_chemistry_is_not_inferred_from_the_label() {
    let missing_type = CIF
        .replace("_atom_site_type_symbol\n", "")
        .replace("C1 C 0.1", "C1 0.1")
        .replace("O1 O 0.2", "O1 0.2");
    assert_eq!(
        parse_small(&missing_type),
        Err(SmallCifError::AtomTypeSymbol)
    );
}

#[test]
fn incomplete_fractional_triplets_are_refused() {
    let broken = CIF
        .replace("_atom_site_fract_z\n", "")
        .replace("C1 C 0.1 0.2 0.3 1", "C1 C 0.1 0.2 1")
        .replace("O1 O 0.2 0.2 0.3 0.5", "O1 O 0.2 0.2 0.5");
    assert_eq!(parse_small(&broken), Err(SmallCifError::Coordinates));
}

#[test]
fn external_cod_fixture_when_configured() {
    let Some(path) = std::env::var_os("MOLFRAME_SMALL_CIF_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).expect("configured COD fixture should be readable");
    let input = InputBuffer::from_bytes(bytes);
    let (document, _) = parse(&input).expect("COD CIF should parse");
    let structure = lower_small_cif_with_options(&document, SmallCifOptions::ddl1())
        .expect("COD CIF should lower");
    assert!(!structure.atoms.is_empty());
    assert!(structure.cell.is_some());
}
