use super::*;
use crate::{Namespace, ReadOptions, read_bytes};

const SOURCE: &str = "data_q\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\nATOM 2 C C2 LIG A 1 1 0 0\n";

fn structure() -> crate::Structure {
    match read_bytes(
        SOURCE.as_bytes().to_vec(),
        Some("query.cif"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

#[test]
fn structure_text_selection_connects_query_and_spatial_execution() {
    let policy = AnalysisPolicy {
        identifiers: Namespace::Label,
        ..AnalysisPolicy::default()
    };
    let selected = structure().select_text("within 1.1 of name C1", &policy);
    let selected = match selected {
        Ok(selected) => selected,
        Err(findings) => panic!("select failed: {findings:?}"),
    };
    assert_eq!(selected.selection.iter().collect::<Vec<_>>(), vec![0, 1]);
}
