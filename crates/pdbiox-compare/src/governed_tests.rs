use super::{
    governed_cad_score, governed_comparison_workflow, governed_contact_map_similarity,
    governed_lddt,
};
use crate::{ContactArea, LddtOptions, PointMapping, PointMatch};
use pdbiox_core::contract::{AlignmentPolicy, AnalysisPolicy};

#[test]
fn governed_score_records_algorithm_parameters_and_coverage() {
    let coordinates = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let options = LddtOptions::standard(15.0);
    let Ok(result) = governed_lddt(
        &coordinates,
        &coordinates,
        &options,
        &AnalysisPolicy::default(),
    ) else {
        panic!("identical coordinates are comparable");
    };
    assert_eq!(result.value.to_bits(), 1.0_f64.to_bits());
    assert_eq!(result.coverage.used, 2);
    assert_eq!(
        result
            .provenance
            .algorithm
            .as_ref()
            .map(pdbiox_core::contract::AlgorithmId::name),
        Some("lddt")
    );
    assert!(result.provenance.parameters.contains_key("tolerances"));
}

#[test]
fn governed_contact_metrics_record_their_input_domain() {
    let contacts = [ContactArea {
        first: 1,
        second: 2,
        area: 4.0,
    }];
    let Ok(cad) = governed_cad_score(&contacts, &contacts, &AnalysisPolicy::default()) else {
        panic!("valid contact areas should score");
    };
    let Ok(similarity) =
        governed_contact_map_similarity(&[(1, 2), (2, 3)], &[(2, 1)], &AnalysisPolicy::default())
    else {
        panic!("contact pairs fit public coverage");
    };
    assert_eq!(cad.value.score.to_bits(), 1.0_f64.to_bits());
    assert_eq!(cad.coverage.used, 1);
    assert_eq!(similarity.value.shared, 1);
    assert_eq!(similarity.coverage.intended, 2);
}

#[test]
fn explicit_alignment_policy_drives_the_workflow_fit() {
    let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let model = [[4.0, 5.0, 0.0], [5.0, 5.0, 0.0], [4.0, 6.0, 0.0]];
    let Ok(mapping) = PointMapping::new(
        (0..3).map(|index| PointMatch {
            reference: index,
            model: index,
        }),
        3,
        3,
    ) else {
        panic!("one-to-one mapping should be valid");
    };
    let policy = AnalysisPolicy {
        alignment: AlignmentPolicy::Explicit("all".into()),
        ..AnalysisPolicy::default()
    };
    let Ok(result) = governed_comparison_workflow(&reference, &model, &mapping, 0.01, &policy)
    else {
        panic!("explicit mapping should fit");
    };
    assert!(result.value.0.rmsd < 1.0e-6);
    assert!(result.value.1.passed);
}
