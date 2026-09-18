use super::{UpdatingSelection, UpdatingSelectionError};
use crate::Timestep;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{ExecutionContext, MemoryBudget, ScratchPolicy};
use molframe_query::{Groups, Query};
use molframe_spatial::SpatialBackend;

fn topology() -> molframe_core::structure::Structure {
    let source = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.auth_atom_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG C1 A 1 0 0 0\nATOM 2 C C2 LIG C2 A 1 1 0 0\n\
ATOM 3 C C3 LIG C3 A 1 4 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    molframe_cif::read(&input, &ReadOptions::new()).map_or_else(
        |findings| panic!("fixture failed: {findings:?}"),
        |(structure, _)| structure,
    )
}

#[test]
fn geometric_membership_updates_without_recompiling_the_query() {
    let query = Query::compile("within 1.1 of label_name C1")
        .unwrap_or_else(|findings| panic!("query failed: {findings:?}"));
    let updating = UpdatingSelection::new(
        topology(),
        &query,
        AnalysisPolicy::default(),
        Groups::new(),
        SpatialBackend::CellList,
    );
    let first = Timestep {
        positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
        ..Timestep::default()
    };
    let second = Timestep {
        positions: vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        ..Timestep::default()
    };
    let selected = |frame| {
        updating
            .evaluate(frame, &ExecutionContext::default())
            .map_or_else(
                |error| panic!("selection failed: {error}"),
                |evaluation| evaluation.selection.iter().collect::<Vec<_>>(),
            )
    };
    assert_eq!(selected(&first), [0, 1]);
    assert_eq!(selected(&second), [0, 2]);
}

#[test]
fn coordinate_rebinding_refuses_before_copying_past_the_shared_budget() {
    let query = Query::compile("all").expect("query");
    let updating = UpdatingSelection::new(
        topology(),
        &query,
        AnalysisPolicy::default(),
        Groups::new(),
        SpatialBackend::CellList,
    );
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(1).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let frame = Timestep {
        positions: vec![[0.0; 3]; 3],
        ..Timestep::default()
    };
    assert!(matches!(
        updating.evaluate(&frame, &context),
        Err(UpdatingSelectionError::Memory(_))
    ));
    assert_eq!(context.reserved_bytes(), 0);
}
