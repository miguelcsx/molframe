use super::*;
use pdbiox_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "data_source\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
_atom_site.pdbx_PDB_model_num\nATOM 44 C CA GLY LONG_CHAIN 9 1 2 3 7\n\
ATOM 45 N N GLY LONG_CHAIN 9 2 2 3 7\n";

#[test]
fn absent_block_identity_refuses_without_an_output_value() {
    let structure = pdbiox_core::Structure::new(pdbiox_core::StructureData::empty());
    assert_eq!(
        write_canonical(&structure),
        Err(CifWriteError::MissingBlockId)
    );
}

#[test]
fn explicit_block_identity_preserves_atom_model_and_chain_identifiers() {
    let structure = read(SOURCE);
    let options = CifWriteOptions::new().with_block_id("chosen");
    let Ok(output) = write_canonical_with_options(&structure, &options) else {
        panic!("complete structure should write");
    };
    assert!(output.starts_with("data_chosen\n"), "{output}");
    assert!(output.contains("ATOM 44 C CA . GLY LONG_CHAIN . 9 ? 1.000 2.000 3.000"));
    assert!(output.ends_with(" 7\n#\n"), "{output}");
}

#[test]
fn absent_atom_site_identity_is_never_replaced_by_an_ordinal() {
    let structure = read(&SOURCE.replace("ATOM 44", "ATOM ?"));
    let options = CifWriteOptions::new().with_block_id("chosen");
    assert_eq!(
        write_canonical_with_options(&structure, &options),
        Err(CifWriteError::MissingAtomSiteId { atom: 0 })
    );
}

#[test]
fn bond_identifiers_require_explicit_generation_permission() {
    let source = format!(
        "{SOURCE}loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
         _struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
         _struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
         _struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
         _struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
         c1 covale LONG_CHAIN 9 GLY CA LONG_CHAIN 9 GLY N\n"
    );
    let structure = read(&source);
    let options = CifWriteOptions::new().with_block_id("chosen");
    assert_eq!(
        write_canonical_with_options(&structure, &options),
        Err(CifWriteError::ConnectionIdsNotEnabled)
    );
}

fn read(source: &str) -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match crate::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}
