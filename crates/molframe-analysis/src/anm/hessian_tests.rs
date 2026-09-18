use super::Hessian;
use crate::network::{ContactGraph, NetworkBudget};
use molframe_core::{ExecutionContext, selection::AtomSelection};
use molframe_spatial::SpatialBackend;

fn budget() -> NetworkBudget {
    NetworkBudget {
        contact_distance: 1.5,
        memory_limit_bytes: 1 << 20,
        backend: SpatialBackend::BruteForce,
    }
}

fn network(positions: &[[f32; 3]]) -> (ContactGraph, Vec<u32>) {
    let Ok(count) = u32::try_from(positions.len()) else {
        panic!("fixture fits an atom index")
    };
    let sites = AtomSelection::from_sorted((0..count).collect());
    let selected: Vec<u32> = (&sites).into_iter().collect();
    let Ok(graph) = ContactGraph::build(
        positions,
        &sites,
        &selected,
        budget(),
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("fixture network builds")
    };
    (graph, selected)
}

#[test]
fn a_single_spring_resists_only_motion_along_its_own_axis() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let (graph, selected) = network(&positions);
    let Ok(hessian) = Hessian::build(&graph, &positions, &selected, budget()) else {
        panic!("fixture Hessian builds")
    };

    // Pulling the two sites apart along the bond stores energy.
    let stretch = [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let mut out = vec![0.0; 6];
    hessian.accumulate(&mut out, &stretch, 1.0);
    assert!(
        (out[0] - (-2.0)).abs() < 1e-12 && (out[3] - 2.0).abs() < 1e-12,
        "stretching gave {out:?}"
    );

    // Sliding them past each other, perpendicular to the bond, does not.
    let shear = [0.0, 1.0, 0.0, 0.0, -1.0, 0.0];
    let mut out = vec![0.0; 6];
    hessian.accumulate(&mut out, &shear, 1.0);
    for value in out {
        assert!(value.abs() < 1e-12, "shear stored energy: {value}");
    }
}

#[test]
fn the_operator_is_symmetric() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.1, 0.0],
        [0.2, 1.0, 0.3],
        [1.1, 1.0, 0.1],
    ];
    let (graph, selected) = network(&positions);
    let Ok(hessian) = Hessian::build(&graph, &positions, &selected, budget()) else {
        panic!("fixture Hessian builds")
    };
    let dimension = hessian.dimension();
    let mut dense = vec![0.0_f64; dimension * dimension];
    let mut unit = vec![0.0_f64; dimension];
    for index in 0..dimension {
        unit.fill(0.0);
        unit[index] = 1.0;
        let mut column = vec![0.0_f64; dimension];
        hessian.accumulate(&mut column, &unit, 1.0);
        for (row, value) in column.iter().enumerate() {
            dense[row * dimension + index] = *value;
        }
    }
    for row in 0..dimension {
        for column in 0..dimension {
            let upper = dense[row * dimension + column];
            let lower = dense[column * dimension + row];
            assert!(
                (upper - lower).abs() < 1e-12,
                "asymmetric at {row},{column}"
            );
        }
    }
}

#[test]
fn a_rigid_translation_stores_no_energy() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    let (graph, selected) = network(&positions);
    let Ok(hessian) = Hessian::build(&graph, &positions, &selected, budget()) else {
        panic!("fixture Hessian builds")
    };
    let field: Vec<f64> = (0..positions.len())
        .flat_map(|_| [0.5_f64, -0.25, 0.75])
        .collect();
    let mut forces = vec![0.0; field.len()];
    hessian.accumulate(&mut forces, &field, 1.0);
    for value in forces {
        assert!(value.abs() < 1e-12, "translation produced force {value}");
    }
}

#[test]
fn coincident_contacts_are_refused_rather_than_normalised() {
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let (graph, selected) = network(&positions);
    assert!(
        Hessian::build(&graph, &positions, &selected, budget()).is_err(),
        "a zero-length spring has no direction to pull along"
    );
}
