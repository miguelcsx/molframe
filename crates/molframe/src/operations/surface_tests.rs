use crate::{
    CoordinateInput, FloatInput, MaskInput, Plan, PlanOperation, PlanValue, SurfaceRequest,
    SurfaceValue,
};

#[test]
fn surface_operations_match_the_facade_surface_kernels() {
    let positions = [[0.0_f32, 0.0, 0.0], [2.5, 0.0, 0.0], [20.0, 0.0, 0.0]];
    let radii = [1.5_f32, 1.5, 1.5];
    let first = [true, false, false];
    let arrays = [CoordinateInput {
        positions: &positions,
    }];
    let floats = [FloatInput { values: &radii }];
    let masks = [MaskInput { values: &first }];
    let probe = 1.4_f32;
    let sample_points = 400_u16;
    let direct_sasa = molframe_surface::shrake_rupley(
        &positions,
        &radii,
        probe,
        sample_points,
        &molframe_core::ExecutionContext::default(),
    )
    .expect("direct SASA");
    let direct_buried = molframe_surface::buried_surface(
        &positions,
        &radii,
        probe,
        sample_points,
        &first,
        &molframe_core::ExecutionContext::default(),
    )
    .expect("direct buried surface");

    let mut plan = Plan::new();
    plan.add(
        "sasa",
        PlanOperation::Surface(Box::new(SurfaceRequest::SolventAccessibleSurface {
            positions: 0,
            radii: 0,
            probe,
            sample_points,
        })),
    )
    .expect("SASA operation");
    plan.add(
        "buried",
        PlanOperation::Surface(Box::new(SurfaceRequest::BuriedSurface {
            positions: 0,
            radii: 0,
            first: 0,
            probe,
            sample_points,
        })),
    )
    .expect("buried surface operation");
    let result = plan
        .execute(
            crate::PlanInput {
                arrays: &arrays,
                floats: &floats,
                masks: &masks,
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native surface plan");
    let values = result
        .entries
        .into_iter()
        .map(|entry| (entry.id.to_string(), entry.value))
        .collect::<std::collections::BTreeMap<_, _>>();
    assert!(matches!(
        values.get("sasa"),
        Some(PlanValue::Surface(value))
            if matches!(value.as_ref(), SurfaceValue::SolventAccessibleSurface(actual) if actual == &direct_sasa)
    ));
    assert!(matches!(
        values.get("buried"),
        Some(PlanValue::Surface(value))
            if matches!(value.as_ref(), SurfaceValue::BuriedSurface(actual) if actual == &direct_buried)
    ));
}
