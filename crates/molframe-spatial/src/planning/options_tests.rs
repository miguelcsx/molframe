use super::*;

#[test]
fn automatic_decisions_are_reproducible_from_the_public_profile() {
    let profile = AutoBackendProfile {
        brute_force_pair_limit: 100,
        kd_target_minimum: 20,
        kd_query_ratio: 4,
        periodic_backend: SpatialBackend::CellList,
    };

    assert_eq!(
        profile.resolve(10, 10, false),
        Ok(SpatialBackend::BruteForce)
    );
    assert_eq!(profile.resolve(4, 40, false), Ok(SpatialBackend::KdTree));
    assert_eq!(profile.resolve(20, 20, false), Ok(SpatialBackend::CellList));
    assert_eq!(profile.resolve(1, 1, true), Ok(SpatialBackend::CellList));
}

#[test]
fn invalid_profiles_name_the_field_that_prevents_a_plan() {
    let options = SpatialSearchOptions {
        automatic: AutoBackendProfile {
            kd_query_ratio: 0,
            ..AutoBackendProfile::BALANCED
        },
        ..SpatialSearchOptions::BALANCED
    };

    assert_eq!(
        options.plan(1, 1, false, 1.0),
        Err(SpatialError::InvalidOption(SpatialOption::KdQueryRatio))
    );
}

#[test]
fn a_neighbor_plan_exposes_its_derived_skin() {
    let options = SpatialSearchOptions {
        backend: SpatialBackend::NeighborList,
        neighbor_skin: NeighborSkinProfile {
            cutoff_ratio: 0.25,
            minimum: 0.75,
        },
        ..SpatialSearchOptions::BALANCED
    };

    let plan = match options.plan(100, 100, false, 4.0) {
        Ok(plan) => plan,
        Err(error) => panic!("plan failed: {error}"),
    };
    assert_eq!(plan.requested_backend, SpatialBackend::NeighborList);
    assert_eq!(plan.backend, SpatialBackend::NeighborList);
    assert_eq!(plan.neighbor_skin, Some(1.0));
}
