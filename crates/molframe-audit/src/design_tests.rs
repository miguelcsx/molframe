use super::*;

fn grid() -> (Vec<PolicyField>, Vec<Vec<usize>>) {
    (
        vec![PolicyField::Altloc, PolicyField::Hydrogens],
        vec![vec![0, 0], vec![0, 1], vec![1, 0], vec![1, 1]],
    )
}

fn scalar(outcomes: [f64; 4]) -> Decomposition {
    let (fields, coordinates) = grid();
    decompose(&fields, &coordinates, |first, second| {
        (outcomes[first] - outcomes[second]).abs()
    })
}

fn close(left: f64, right: f64) -> bool {
    (left - right).abs() < 1e-12
}

#[test]
fn additive_decisions_leave_nothing_to_the_interaction() {
    let result = scalar([0.0, 1.0, 10.0, 11.0]);
    assert!(close(result.total_variation, 101.0));
    assert!(close(result.main_effects[0].share, 100.0 / 101.0));
    assert!(close(result.main_effects[1].share, 1.0 / 101.0));
    assert!(close(result.interactions[0].share, 0.0));
    assert!(close(result.higher_order, 0.0));
    assert_eq!(result.interactions[0].first, PolicyField::Altloc);
    assert_eq!(result.interactions[0].second, PolicyField::Hydrogens);
}

#[test]
fn a_decision_that_matters_only_with_another_shows_up_as_an_interaction() {
    // Only the combination (1, 1) differs: each decision alone explains a third.
    let result = scalar([0.0, 0.0, 0.0, 10.0]);
    assert!(close(result.total_variation, 75.0));
    assert!(close(result.main_effects[0].share, 1.0 / 3.0));
    assert!(close(result.main_effects[1].share, 1.0 / 3.0));
    assert!(close(result.interactions[0].share, 1.0 / 3.0));
    // Changing one decision with the other held at its first level moves nothing.
    assert!(close(result.main_effects[0].mean_change, 5.0));
}

#[test]
fn identical_runs_explain_nothing() {
    let result = scalar([2.0; 4]);
    assert!(close(result.total_variation, 0.0));
    assert!(close(result.mean_distance, 0.0));
    assert!(close(result.max_distance, 0.0));
    assert!(
        result
            .main_effects
            .iter()
            .all(|effect| close(effect.share, 0.0))
    );
    assert!(close(result.higher_order, 0.0));
}

#[test]
fn three_way_structure_is_reported_as_higher_order() {
    let fields = vec![
        PolicyField::Altloc,
        PolicyField::Hydrogens,
        PolicyField::Model,
    ];
    let coordinates: Vec<Vec<usize>> = (0..8)
        .map(|run| vec![(run >> 2) & 1, (run >> 1) & 1, run & 1])
        .collect();
    // The parity of the three levels is zero in every main effect and pair.
    let outcome = |run: usize| -> f64 {
        if run.count_ones().is_multiple_of(2) {
            0.0
        } else {
            1.0
        }
    };
    let result = decompose(&fields, &coordinates, |first, second| {
        (outcome(first) - outcome(second)).abs()
    });
    assert!(
        result
            .main_effects
            .iter()
            .all(|effect| close(effect.share, 0.0))
    );
    assert!(
        result
            .interactions
            .iter()
            .all(|effect| close(effect.share, 0.0))
    );
    assert!(close(result.higher_order, 1.0));
}

fn share_of(result: &Decomposition, field: PolicyField) -> f64 {
    result
        .shapley
        .iter()
        .find(|attribution| attribution.key == field)
        .map_or(f64::NAN, |attribution| attribution.share)
}

#[test]
fn for_additive_decisions_the_shapley_share_is_the_main_effect() {
    let result = scalar([0.0, 1.0, 10.0, 11.0]);
    assert!(result.balanced);
    assert!(close(share_of(&result, PolicyField::Altloc), 100.0 / 101.0));
    assert!(close(
        share_of(&result, PolicyField::Hydrogens),
        1.0 / 101.0
    ));
}

#[test]
fn a_pure_interaction_is_shared_equally_between_the_decisions_that_make_it() {
    let result = scalar([0.0, 0.0, 0.0, 10.0]);
    assert!(close(share_of(&result, PolicyField::Altloc), 0.5));
    assert!(close(share_of(&result, PolicyField::Hydrogens), 0.5));
    let total: f64 = result
        .shapley
        .iter()
        .map(|attribution| attribution.share)
        .sum();
    assert!(close(total, 1.0));
}

#[test]
fn shapley_shares_exist_and_sum_to_one_when_a_combination_is_forbidden() {
    // The (1, 1) cell is missing: the design is not the whole product, so the additive
    // main-effect and interaction shares do not apply, and the Shapley shares still do.
    let fields = vec![PolicyField::Altloc, PolicyField::Hydrogens];
    let coordinates = vec![vec![0, 0], vec![0, 1], vec![1, 0]];
    let outcomes = [0.0_f64, 4.0, 1.0];
    let result = decompose(&fields, &coordinates, |first, second| {
        (outcomes[first] - outcomes[second]).abs()
    });
    assert!(!result.balanced);
    let total: f64 = result
        .shapley
        .iter()
        .map(|attribution| attribution.share)
        .sum();
    assert!(close(total, 1.0));
    // Under the uniform distribution on the three valid cells, Var(Y) = 26/9.
    // Holding altloc explains 1/13; holding hydrogens explains 49/52.
    // Averaging the two orders gives 7/104 and 97/104, independently of distances.
    assert!(close(share_of(&result, PolicyField::Altloc), 7.0 / 104.0));
    assert!(close(
        share_of(&result, PolicyField::Hydrogens),
        97.0 / 104.0
    ));
}

#[test]
fn shares_are_also_given_per_class_of_uncertainty() {
    let fields = vec![
        PolicyField::Model,
        PolicyField::Hydrogens,
        PolicyField::Precision,
    ];
    let coordinates: Vec<Vec<usize>> = (0..8)
        .map(|run| vec![(run >> 2) & 1, (run >> 1) & 1, run & 1])
        .collect();
    // Only the hydrogens decision changes the answer.
    let outcome = |run: usize| -> f64 { if (run >> 1) & 1 == 1 { 5.0 } else { 0.0 } };
    let result = decompose(&fields, &coordinates, |first, second| {
        (outcome(first) - outcome(second)).abs()
    });
    let class = |wanted: UncertaintyClass| {
        result
            .by_class
            .iter()
            .find(|attribution| attribution.key == wanted)
            .map_or(f64::NAN, |attribution| attribution.share)
    };
    assert!(close(class(UncertaintyClass::Interpretive), 1.0));
    assert!(close(class(UncertaintyClass::Structural), 0.0));
    assert!(close(class(UncertaintyClass::Numerical), 0.0));
    assert_eq!(
        UncertaintyClass::of(PolicyField::VdwRadii),
        UncertaintyClass::Algorithmic
    );
}
