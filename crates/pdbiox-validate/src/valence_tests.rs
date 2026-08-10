use super::overvalent_atoms;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

const ATOMS: &str = "data_v\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 O O1 LIG A 1 0 0 0\n\
ATOM 2 C C1 LIG A 1 1 0 0\n\
ATOM 3 C C2 LIG A 1 0 1 0\n\
ATOM 4 C C3 LIG A 1 0 0 1\n";

const CONN_HEADER: &str = "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_value_order\n";

#[test]
fn an_oxygen_with_three_bonds_is_over_coordinated() {
    let source = format!(
        "{ATOMS}{CONN_HEADER}\
1 covale A 1 LIG O1 A 1 LIG C1 SING\n\
2 covale A 1 LIG O1 A 1 LIG C2 SING\n\
3 covale A 1 LIG O1 A 1 LIG C3 SING\n"
    );
    let errors = overvalent_atoms(&structure(&source));
    assert_eq!(errors.len(), 1);
    assert_eq!(errors[0].atom.get(), 0);
    assert_eq!((errors[0].bonds, errors[0].maximum), (3, 2));
}

#[test]
fn an_oxygen_with_two_bonds_is_fine() {
    let source = format!(
        "{ATOMS}{CONN_HEADER}\
1 covale A 1 LIG O1 A 1 LIG C1 SING\n\
2 covale A 1 LIG O1 A 1 LIG C2 SING\n"
    );
    assert!(overvalent_atoms(&structure(&source)).is_empty());
}

#[test]
fn without_bonds_nothing_is_checked() {
    assert!(overvalent_atoms(&structure(ATOMS)).is_empty());
}
