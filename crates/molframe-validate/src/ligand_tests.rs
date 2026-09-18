use super::{ligand_geometry, ligand_geometry_outliers};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

const CONN: &str = "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_value_order\n\
1 covale A 1 LIG C1 A 1 LIG C2 SING\n";

fn ligand(kind: &str, separation: f32) -> String {
    format!(
        "data_l\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
{kind} 1 C C1 LIG A 1 0 0 0\n\
{kind} 2 C C2 LIG A 1 {separation} 0 0\n\
{CONN}"
    )
}

#[test]
fn a_stretched_ligand_bond_is_flagged() {
    let report = ligand_geometry(&structure(&ligand("HETATM", 3.0)), 0.3);
    assert_eq!(report.outliers.len(), 1);
    assert!(report.outliers[0].deviation > 0.3);
    assert_eq!((report.intended, report.assessed), (1, 1));
}

#[test]
fn a_normal_ligand_bond_is_not_flagged() {
    let deviations = ligand_geometry_outliers(&structure(&ligand("HETATM", 1.54)), 0.3);
    assert!(deviations.is_empty());
}

#[test]
fn a_stretched_polymer_bond_is_left_to_the_general_check() {
    // ATOM records are polymer residues, not het, so this ligand check ignores
    // them even when the bond is stretched.
    let deviations = ligand_geometry_outliers(&structure(&ligand("ATOM", 3.0)), 0.3);
    assert!(deviations.is_empty());
}
