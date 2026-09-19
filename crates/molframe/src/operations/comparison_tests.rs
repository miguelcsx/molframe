//! Comparison operations over the plan, on borrowed coordinate inputs.

use super::{
    ComparisonMetric, ComparisonRequest, CoordinateInput, CoordinateSlot, Plan, PlanValue,
};
use crate::PlanInput;

#[test]
fn comparison_operations_share_borrowed_coordinate_inputs() {
    let mobile = [[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let reference = mobile;
    let arrays = [
        CoordinateInput { positions: &mobile },
        CoordinateInput {
            positions: &reference,
        },
    ];
    let mut plan = Plan::new();
    plan.add(
        "lddt",
        ComparisonRequest::new(
            CoordinateSlot::new(0),
            CoordinateSlot::new(1),
            ComparisonMetric::Lddt {
                inclusion_radius: 15.0,
            },
        ),
    )
    .expect("lDDT operation");
    plan.add(
        "tm",
        ComparisonRequest::new(
            CoordinateSlot::new(0),
            CoordinateSlot::new(1),
            ComparisonMetric::TmScore,
        ),
    )
    .expect("TM-score operation");
    plan.add(
        "gdt",
        ComparisonRequest::new(
            CoordinateSlot::new(0),
            CoordinateSlot::new(1),
            ComparisonMetric::GdtTs,
        ),
    )
    .expect("GDT operation");

    let result = plan
        .execute(
            PlanInput {
                structure: None,
                arrays: &arrays,
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native comparison plan");
    assert_eq!(result.entries.len(), 3);
    assert!(result.entries.iter().all(|entry| {
        matches!(
            entry.value,
            PlanValue::Comparison(value) if (value.value - 1.0).abs() < f64::EPSILON
        )
    }));
}
