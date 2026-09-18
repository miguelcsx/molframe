use super::{CrystalNeighborOptions, collect_crystal_neighbors};
use crate::lower_symmetry;
use molframe_cif::{parse, read};
use molframe_core::{Code, ExecutionContext, InputBuffer, ModelIndex, ReadOptions, Structure};
use molframe_spatial::SpatialBackend;

const ATOM_AND_CELL: &str = r"data_crystal
_cell.length_a 10
_cell.length_b 100
_cell.length_c 100
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 1 A 1 0 0
loop_
_space_group_symop.id
_space_group_symop.operation_xyz
1 'x,y,z'
";

#[test]
fn identity_image_is_excluded_but_one_periodic_pair_is_kept() {
    let (structure, symmetry) = fixture(ATOM_AND_CELL);
    let neighbors = match collect_crystal_neighbors(
        &structure,
        &symmetry,
        ModelIndex::new(0),
        10.1,
        CrystalNeighborOptions::default(),
        &ExecutionContext::default(),
    ) {
        Ok(neighbors) => neighbors,
        Err(finding) => panic!("search failed: {finding}"),
    };
    assert_eq!(neighbors.len(), 1);
    assert_eq!(neighbors[0].lattice, [-1, 0, 0]);
    assert!((neighbors[0].distance_squared - 100.0).abs() < 1e-9);
    for backend in [
        SpatialBackend::BruteForce,
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
    ] {
        let actual = collect_crystal_neighbors(
            &structure,
            &symmetry,
            ModelIndex::new(0),
            10.1,
            CrystalNeighborOptions {
                backend,
                ..CrystalNeighborOptions::default()
            },
            &ExecutionContext::default(),
        );
        assert_eq!(actual.ok().as_deref(), Some(neighbors.as_slice()));
    }
}

#[test]
fn self_inverse_operator_is_a_contact_when_positions_differ() {
    let text = ATOM_AND_CELL.replace("1 'x,y,z'", "1 'x,y,z'\n2 '-x,-y,-z'");
    let (structure, symmetry) = fixture(&text);
    let neighbors = match collect_crystal_neighbors(
        &structure,
        &symmetry,
        ModelIndex::new(0),
        2.1,
        CrystalNeighborOptions::default(),
        &ExecutionContext::default(),
    ) {
        Ok(neighbors) => neighbors,
        Err(finding) => panic!("search failed: {finding}"),
    };
    assert_eq!(neighbors.len(), 1);
    assert_eq!(neighbors[0].operation, 1);
    assert!((neighbors[0].distance_squared - 4.0).abs() < 1e-9);
}

#[test]
fn invalid_inputs_and_candidate_ceiling_are_diagnostics() {
    let (structure, symmetry) = fixture(ATOM_AND_CELL);
    let invalid = collect_crystal_neighbors(
        &structure,
        &symmetry,
        ModelIndex::new(0),
        0.0,
        CrystalNeighborOptions::default(),
        &ExecutionContext::default(),
    );
    assert_eq!(
        invalid.err().map(|finding| finding.code()),
        Some(Code::E6016)
    );
    let limited = collect_crystal_neighbors(
        &structure,
        &symmetry,
        ModelIndex::new(0),
        1.0,
        CrystalNeighborOptions {
            candidate_limit: 0,
            ..CrystalNeighborOptions::default()
        },
        &ExecutionContext::default(),
    );
    assert_eq!(
        limited.err().map(|finding| finding.code()),
        Some(Code::E6017)
    );

    let (without_cell, symmetry) = fixture(
        &ATOM_AND_CELL
            .lines()
            .filter(|line| !line.starts_with("_cell."))
            .collect::<Vec<_>>()
            .join("\n"),
    );
    let missing = collect_crystal_neighbors(
        &without_cell,
        &symmetry,
        ModelIndex::new(0),
        1.0,
        CrystalNeighborOptions::default(),
        &ExecutionContext::default(),
    );
    assert_eq!(
        missing.err().map(|finding| finding.code()),
        Some(Code::E5004)
    );
}

fn fixture(text: &str) -> (Structure, crate::SymmetrySet) {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    let symmetry = match lower_symmetry(&document) {
        Ok(symmetry) => symmetry,
        Err(findings) => panic!("symmetry lowering failed: {findings:?}"),
    };
    let structure = match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("structure lowering failed: {findings:?}"),
    };
    (structure, symmetry)
}
