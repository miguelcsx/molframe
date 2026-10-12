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

#[test]
fn facade_selection_uses_the_callers_cancellation_for_simple_and_spatial_queries() {
    let structure = structure();
    let policy = AnalysisPolicy {
        identifiers: Namespace::Label,
        ..AnalysisPolicy::default()
    };
    let context = ExecutionContext::default();
    for query in ["name C1", "within 1.1 of name C1"] {
        let expected = structure.select(query, &policy).expect("default query");
        let actual = structure
            .select_with_context(query, &policy, &context)
            .expect("governed query");
        assert_eq!(
            actual
                .atoms()
                .map(molframe_core::structure::AtomRef::index)
                .collect::<Vec<_>>(),
            expected
                .atoms()
                .map(molframe_core::structure::AtomRef::index)
                .collect::<Vec<_>>()
        );
    }
    context.cancellation().cancel();
    for query in ["name C1", "within 1.1 of name C1"] {
        let failure = structure
            .select_with_context(query, &policy, &context)
            .expect_err("cancelled");
        assert!(failure.to_string().contains("E1904"));
    }
}

const SHELL: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA ALA A 1 0 0 0\nATOM 2 C CB ALA A 1 1 0 0\n\
HETATM 3 C C1 LIG B 2 0 3 0\nATOM 4 C CA GLY A 3 20 0 0\n";

#[test]
fn context_reading_operands_keep_their_written_order_under_cost_planning() {
    let structure = match read_bytes(
        SHELL.as_bytes().to_vec(),
        Some("shell.cif"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    let policy = AnalysisPolicy {
        identifiers: Namespace::Label,
        ..AnalysisPolicy::default()
    };
    let rows = |text: &str| match structure.select_text(text, &policy) {
        Ok(selected) => selected.selection.iter().collect::<Vec<_>>(),
        Err(findings) => panic!("select failed: {findings:?}"),
    };
    // The membership test is cheaper than the geometric one, so a cost-ordering
    // planner would run it first and hide the ligand from `within`.
    let shell = rows("(within 4 of (resname LIG)) and (not (resname LIG))");
    assert_eq!(shell, vec![0, 1]);
    assert_eq!(
        rows("(not (resname LIG)) and (within 4 of (resname LIG))"),
        Vec::<u32>::new(),
        "the written order is the evaluation order for context-reading operands"
    );
}

const TWO_CHAINS: &str = "data_q\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A H 1 0 0 0\nATOM 2 C C2 LIG A H 1 1 0 0\n\
ATOM 3 C C3 LIG B L 1 2 0 0\nATOM 4 C C4 LIG B L 1 50 0 0\n";

#[test]
fn a_conjunction_selects_the_same_atoms_in_either_order() {
    let structure = match read_bytes(
        TWO_CHAINS.as_bytes().to_vec(),
        Some("two.cif"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    let policy = AnalysisPolicy::default();
    let pick = |query: &str| -> Vec<u32> {
        match structure.select_text(query, &policy) {
            Ok(selected) => selected.selection.iter().collect(),
            Err(findings) => panic!("select failed: {findings:?}"),
        }
    };
    let forward = pick("chain H and within 1.5 of chain L");
    assert_eq!(forward, vec![1]);
    assert_eq!(forward, pick("within 1.5 of chain L and chain H"));
    assert_eq!(
        pick("chain H and byres chain L"),
        pick("byres chain L and chain H")
    );
}
