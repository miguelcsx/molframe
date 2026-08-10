use super::bond_length_deviations;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
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

fn two_carbons(separation: f32) -> String {
    format!(
        "data_b\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 {separation} 0 0\n\
{CONN}"
    )
}

#[test]
fn a_normal_carbon_carbon_bond_is_not_flagged() {
    // ~1.54 Å is a standard C-C bond; the covalent-radii sum is ~1.5 Å.
    let deviations = bond_length_deviations(&structure(&two_carbons(1.54)), 0.3);
    assert!(deviations.is_empty(), "normal bond flagged: {deviations:?}");
}

#[test]
fn a_stretched_bond_is_flagged() {
    let deviations = bond_length_deviations(&structure(&two_carbons(3.0)), 0.3);
    assert_eq!(deviations.len(), 1);
    assert_eq!(
        (deviations[0].atom_a.get(), deviations[0].atom_b.get()),
        (0, 1)
    );
    assert!(deviations[0].deviation > 0.3, "expected a positive stretch");
}

#[test]
fn without_bonds_nothing_is_measured() {
    let no_conn = "data_b\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 3 0 0\n";
    assert!(bond_length_deviations(&structure(no_conn), 0.3).is_empty());
}
