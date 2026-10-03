use super::*;
use molframe_core::io::{InputBuffer, ReadOptions};

fn structure() -> Structure {
    let text = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C C1 LIG A 1 0 0 0\n";
    molframe_cif::read(
        &InputBuffer::from_bytes(text.as_bytes().to_vec()),
        &ReadOptions::new(),
    )
    .unwrap()
    .0
}

fn spec(origin: f64, spacing: f64, count: usize) -> GridSpec {
    GridSpec {
        dimensions: [count, 1, 1],
        voxel_to_world: [
            [spacing, 0.0, 0.0, origin],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    }
}

#[test]
fn a_unit_charge_matches_the_thermal_coulomb_definition_at_five_and_ten_angstrom() {
    let field = contact_potential(&structure(), &[1.0], spec(5.0, 5.0, 2), 12.0).unwrap();
    // Independent SI definition: potential in volts divided by kT/e.
    let elementary = 1.602_176_634e-19;
    let voltage_unit = 1.380_649e-23 * 298.0 / elementary;
    for (value, radius) in field.values.iter().zip([5.0_f64, 10.0]) {
        let volts = elementary
            / (4.0
                * std::f64::consts::PI
                * 8.854_187_812_8e-12
                * (4.0 * radius)
                * radius
                * 1.0e-10);
        let expected = volts / voltage_unit;
        assert!((value / expected - 1.0).abs() < 1.0e-4);
    }
}

#[test]
fn atom_centres_are_softened_and_the_finite_cutoff_is_inclusive() {
    let field = contact_potential(&structure(), &[-1.0], spec(0.0, 0.5, 3), 12.0).unwrap();
    assert!(field.values[0].is_finite() && field.values[0] < 0.0);
    assert_eq!(field.values[0].to_bits(), field.values[1].to_bits());
    assert_eq!(field.values[0].to_bits(), field.values[2].to_bits());
    let edge = contact_potential(&structure(), &[1.0], spec(12.0, 0.01, 2), 12.0).unwrap();
    assert!(edge.values[0] > 0.0);
    assert_eq!(edge.values[1].to_bits(), 0.0_f64.to_bits());
}

#[test]
fn skew_affines_and_storage_order_match_direct_distances() {
    let mut grid = spec(2.0, 1.0, 2);
    grid.dimensions = [2, 2, 2];
    grid.voxel_to_world[0][1] = 0.5;
    let field = contact_potential(&structure(), &[1.0], grid, 12.0).unwrap();
    for z in 0..2 {
        for y in 0..2 {
            for x in 0..2 {
                let radius_squared =
                    (2.0 + f64::from(x) + 0.5 * f64::from(y)).powi(2) + f64::from(y * y + z * z);
                let index = usize::try_from((z * 2 + y) * 2 + x).unwrap();
                assert!(
                    (field.values[index] - THERMAL_COULOMB / (4.0 * radius_squared)).abs()
                        < 1.0e-10
                );
            }
        }
    }
}

#[test]
fn malformed_geometry_and_charge_columns_are_rejected() {
    assert_eq!(
        contact_potential(&structure(), &[], spec(0.0, 1.0, 1), 12.0),
        Err(PotentialError::MisalignedInput)
    );
    assert_eq!(
        contact_potential(&structure(), &[f64::NAN], spec(0.0, 1.0, 1), 12.0),
        Err(PotentialError::NonFiniteInput)
    );
    assert_eq!(
        contact_potential(&structure(), &[1.0], spec(0.0, 0.0, 1), 12.0),
        Err(PotentialError::InvalidGrid)
    );
    assert_eq!(
        contact_potential(&structure(), &[1.0], spec(0.0, 1.0, 0), 12.0),
        Err(PotentialError::InvalidGrid)
    );
}

#[test]
fn fixed_voxel_blocks_have_identical_values_at_every_worker_budget() {
    let structure = structure();
    let grid = spec(-20.0, 0.001, 16_385);
    let serial = contact_potential_in(
        &structure,
        &[1.0],
        grid,
        12.0,
        &ExecutionContext::builder()
            .worker_budget(1)
            .build()
            .unwrap(),
    )
    .unwrap();
    for workers in [2, 4] {
        let parallel = contact_potential_in(
            &structure,
            &[1.0],
            grid,
            12.0,
            &ExecutionContext::builder()
                .worker_budget(workers)
                .build()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(serial, parallel);
    }
}

#[test]
fn affine_points_near_the_cutoff_are_not_lost_to_f32_candidate_rounding() {
    let field = contact_potential(
        &structure(),
        &[1.0],
        spec(11.999_999_9, 0.000_000_2, 2),
        12.0,
    )
    .unwrap();
    assert!(field.values[0] > 0.0);
    assert_eq!(field.values[1].to_bits(), 0.0_f64.to_bits());
}
