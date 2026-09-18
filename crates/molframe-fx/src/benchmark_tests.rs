use super::{AmeMeasurementError, measure_ame, measure_motifbench};

const SHAPE: [[f32; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn shifted(offset: [f32; 3]) -> Vec<[f32; 3]> {
    SHAPE
        .iter()
        .map(|point| {
            [
                point[0] + offset[0],
                point[1] + offset[1],
                point[2] + offset[2],
            ]
        })
        .collect()
}

#[test]
fn motifbench_fits_scaffold_and_motif_independently() {
    let scaffold = shifted([4.0, -2.0, 1.0]);
    let motif = shifted([-3.0, 5.0, 2.0]);
    let measured = measure_motifbench(&scaffold, &SHAPE, &SHAPE, &motif)
        .expect("non-degenerate coordinates should fit");
    assert!(measured.rmsd < 1e-6);
    assert!(measured.motif_rmsd < 1e-6);
}

#[test]
fn ame_uses_backbone_fit_for_heavy_atoms_and_prediction_frame_for_clashes() {
    let generated_backbone = shifted([3.0, 2.0, -1.0]);
    let generated_heavy = shifted([3.0, 2.0, -1.0]);
    let measured = measure_ame(
        &generated_backbone,
        &SHAPE,
        &generated_heavy,
        &SHAPE,
        &[[0.0, 0.0, 2.0]],
        &[[0.0, 0.0, 0.0], [8.0, 0.0, 0.0]],
    )
    .expect("valid AME atom sets should measure");
    assert!(measured.catalytic_heavy_atom_rmsd < 1e-6);
    assert!((measured.ligand_backbone_min_distance - 2.0).abs() < 1e-6);
}

#[test]
fn ame_refuses_an_empty_ligand_instead_of_inventing_a_distance() {
    let error = measure_ame(&SHAPE, &SHAPE, &SHAPE, &SHAPE, &[], &SHAPE)
        .expect_err("an empty ligand has no clash distance");
    assert_eq!(error, AmeMeasurementError::EmptyLigand);
}
