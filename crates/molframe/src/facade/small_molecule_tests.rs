//! SDF, MOL2 and core-CIF reads and the SDF write through the facade.

use super::*;

const WATER_SDF: &str = "water\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n\
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0\n  1  3  1  0\nM  END\n$$$$\n";

const WATER_MOL2: &str = "@<TRIPOS>MOLECULE\nwat\n 3 2 1 0 0\nSMALL\nNO_CHARGES\n\n\n\
@<TRIPOS>ATOM\n\
 1 OW 0.0 0.0 0.0 O.3 1 WAT 0.0\n\
 2 HW1 0.757 0.586 0.0 H 1 WAT 0.0\n\
 3 HW2 -0.757 0.586 0.0 H 1 WAT 0.0\n\
@<TRIPOS>BOND\n 1 1 2 1\n 2 1 3 1\n";

const SALT: &str = "data_salt\n_cell_length_a 5.64\n_cell_length_b 5.64\n_cell_length_c 5.64\n\
_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
Na1 Na 0.0 0.0 0.0\nCl1 Cl 0.5 0.5 0.5\n";

fn read_text(text: &str, name: &str) -> (Structure, Vec<Diagnostic>) {
    read_bytes(text.as_bytes().to_vec(), Some(name), &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("{name} failed to read: {findings:?}"))
}

#[test]
fn sdf_and_mol2_read_with_their_atoms_and_bonds() {
    for (text, name) in [(WATER_SDF, "water.sdf"), (WATER_MOL2, "water.mol2")] {
        let (structure, _) = read_text(text, name);
        assert_eq!(structure.atom_count(), 3, "{name}");
        assert_eq!(structure.engine().data().bonds.len(), 2, "{name}");
    }
}

#[test]
fn small_molecule_content_is_detected_without_a_suffix() {
    for text in [WATER_SDF, WATER_MOL2, SALT] {
        let (structure, _) = read_bytes(text.as_bytes().to_vec(), None, &ReadOptions::new())
            .expect("content alone names the format");
        assert!(structure.atom_count() >= 2);
    }
}

#[test]
fn a_core_cif_places_fractional_sites_in_the_cell() {
    let (structure, _) = read_text(SALT, "salt.cif");
    let positions: Vec<_> = structure
        .engine()
        .data()
        .atoms()
        .filter_map(molframe_core::structure::AtomRef::position)
        .collect();
    assert_eq!(positions, [[0.0; 3], [2.82; 3]]);
}

#[test]
fn a_multi_record_sdf_reads_the_first_and_warns() {
    let two = format!("{WATER_SDF}{WATER_SDF}");
    let (structure, findings) = read_text(&two, "two.sdf");
    assert_eq!(structure.atom_count(), 3);
    assert!(findings.iter().any(|finding| finding.code() == Code::W1001));
}

#[test]
fn an_sdf_write_reads_back_within_coordinate_precision() {
    let (structure, _) = read_text(WATER_SDF, "water.sdf");
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("round.sdf");
    write(&path, &structure).expect("SDF writes");
    let (again, _) = read_with_diagnostics(&path).expect("SDF reads back");
    let positions = |structure: &Structure| -> Vec<[f32; 3]> {
        structure
            .engine()
            .data()
            .atoms()
            .filter_map(molframe_core::structure::AtomRef::position)
            .collect()
    };
    for (before, after) in positions(&structure).iter().zip(positions(&again)) {
        for axis in 0..3 {
            assert!((before[axis] - after[axis]).abs() < 1e-3);
        }
    }
    assert_eq!(again.engine().data().bonds.len(), 2);
}
