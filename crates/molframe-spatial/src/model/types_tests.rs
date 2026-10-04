use super::*;

#[test]
fn resource_and_cancellation_failures_do_not_masquerade_as_bad_selection_values() {
    let memory = SpatialError::Memory(molframe_core::MemoryBudgetError::Exhausted {
        requested: 100,
        available: 10,
    });
    assert_eq!(memory.into_diagnostic().code(), Code::E7001);
    assert_eq!(
        SpatialError::Cancelled.into_diagnostic().code(),
        Code::E1904
    );
}

#[test]
fn every_backend_round_trips_through_its_name() {
    for name in SpatialBackend::NAMES {
        let backend: SpatialBackend = name.parse().expect("a listed name parses");
        assert_eq!(backend.name(), *name);
        assert_eq!(backend.to_string(), *name);
    }
}

#[test]
fn backend_spellings_accept_the_short_and_underscore_forms() {
    assert_eq!("cell".parse(), Ok(SpatialBackend::CellList));
    assert_eq!("kd_tree".parse(), Ok(SpatialBackend::KdTree));
    assert_eq!("brute_force".parse(), Ok(SpatialBackend::BruteForce));
    let refused = "octree"
        .parse::<SpatialBackend>()
        .expect_err("not a backend");
    assert_eq!(refused.field, "backend");
}
