use super::*;
use crate::{Groups, Query};
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;

const SOURCE: &str = "data_plan\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 O O  GLY A 1 1 0 0\n";

#[test]
fn binding_resolves_globs_to_symbols_and_matches_direct_evaluation() {
    let structure = structure();
    let policy =
        AnalysisPolicy::default().with_identifiers(pdbiox_core::contract::Namespace::Label);
    let query = match Query::compile("name C*") {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    let plan = query.plan(&structure, &policy);
    assert!(matches!(plan.expr, PhysicalExpr::ResolvedMembership { .. }));
    let result = match plan.evaluate(&structure, &policy, &Groups::new(), None) {
        Ok(result) => result,
        Err(findings) => panic!("physical evaluation failed: {findings:?}"),
    };
    assert_eq!(result.selection.iter().collect::<Vec<_>>(), vec![0]);
}

#[test]
fn constant_branches_are_folded_before_execution() {
    let structure = structure();
    let policy = AnalysisPolicy::default();
    let query = match Query::compile("all and none") {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    assert!(matches!(
        query.plan(&structure, &policy).expr,
        PhysicalExpr::None
    ));
}

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}
