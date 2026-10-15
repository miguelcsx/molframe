//! Pruned evaluation equals the unpruned evaluation over the admitted mappings.

use super::{residue_motif, seeded};
use crate::{
    Admission, Comparison, Evaluation, MappingPlan, MeasurementOptions, MissingVerdict,
    VerdictProfile, VerdictRule, VerdictStatus, evaluate_plan,
};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::execution::CancellationToken;
use molframe_core::structure::Structure;
use proptest::prelude::*;
use std::sync::Arc;

fn evaluated(
    structure: &Structure,
    motif: &crate::Motif,
    profile: &VerdictProfile,
    admission: Admission<'_>,
) -> Vec<Evaluation> {
    let plan = MappingPlan::build_admitting(
        structure,
        motif,
        None,
        &AnalysisPolicy::default(),
        admission,
    )
    .unwrap_or_else(|error| panic!("plan failed: {error}"));
    evaluate_plan(
        structure,
        motif,
        &Arc::new(plan),
        profile,
        usize::MAX,
        MeasurementOptions {
            maximum_alternatives: 256,
            plane_fit: molframe_geom::EigenOptions::standard(),
        },
        &CancellationToken::new(),
    )
    .unwrap_or_else(|error| panic!("evaluation failed: {error}"))
    .evaluations
}

fn comparison(kind: u8, low: f64, high: f64) -> Comparison {
    match kind {
        0 => Comparison::Between(low, high),
        1 => Comparison::AtMost(high),
        2 => Comparison::AtLeast(low),
        3 => Comparison::LessThan(high),
        _ => Comparison::GreaterThan(low),
    }
}

/// Whether a rule's closed range holds `value`.
fn within(comparison: Comparison, value: f64) -> bool {
    match comparison {
        Comparison::Between(low, high) => value >= low && value <= high,
        Comparison::AtMost(limit) | Comparison::LessThan(limit) => value <= limit,
        Comparison::AtLeast(limit) | Comparison::GreaterThan(limit) => value >= limit,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn pruned_evaluation_is_the_admitted_slice_of_the_unpruned_one(
        seed in 1_u64..u64::MAX,
        targets in (2.0_f64..14.0, 2.0_f64..14.0),
        tolerances in (0.0_f64..3.0, 0.0_f64..3.0),
        windows in (0.0_f64..10.0, 0.0_f64..8.0),
        kinds in (0_u8..5, 0_u8..5),
        ruled in (any::<bool>(), any::<bool>()),
    ) {
        let structure = seeded(45, seed);
        let motif = residue_motif(
            &[("a", "SER", "OG"), ("b", "TRP", "NE1"), ("c", "MET", "SD")],
            &[
                ("a", "b", targets.0, tolerances.0),
                ("b", "c", targets.1, tolerances.1),
            ],
        );
        let (low, high) = (windows.0, windows.0 + windows.1);
        let rules: Vec<_> = [(ruled.0, "d0", kinds.0), (ruled.1, "d1", kinds.1)]
            .into_iter()
            .filter(|(present, ..)| *present)
            .map(|(_, metric, kind)| VerdictRule {
                metric: metric.into(),
                comparison: comparison(kind, low, high),
            })
            .collect();
        let profile = VerdictProfile::new("property-1.0", rules.clone(), MissingVerdict::Fail);
        let full = evaluated(&structure, &motif, &profile, Admission::All);
        let pruned = evaluated(&structure, &motif, &profile, Admission::Passable(&profile));

        let admitted = |evaluation: &Evaluation| {
            evaluation.measurements.constraints.iter().all(|constraint| {
                let value = constraint.value.map(crate::MeasurementValue::numeric);
                let by_rule = value.is_some_and(|value| {
                    rules
                        .iter()
                        .filter(|rule| rule.metric == constraint.name)
                        .all(|rule| within(rule.comparison, value))
                });
                constraint.satisfied == Some(true) || by_rule
            })
        };
        let expected: Vec<_> = full.iter().filter(|evaluation| admitted(evaluation)).collect();
        prop_assert_eq!(pruned.len(), expected.len());
        for (kept, reference) in pruned.iter().zip(expected) {
            prop_assert_eq!(&kept.residues, &reference.residues);
            prop_assert_eq!(&kept.measurements, &reference.measurements);
            prop_assert_eq!(&kept.verdict, &reference.verdict);
        }
        let passing = |evaluations: &[Evaluation]| -> Vec<Box<[u32]>> {
            evaluations
                .iter()
                .filter(|evaluation| evaluation.verdict.status == VerdictStatus::Pass)
                .map(|evaluation| evaluation.residues.clone())
                .collect()
        };
        prop_assert_eq!(passing(&pruned), passing(&full));
    }
}
