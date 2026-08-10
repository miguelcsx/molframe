use super::{
    CandidateAggregation, ProfileAlignment, ProfileAtomSet, ame_heavy_atom_1_0, motifbench_1_0,
};
use crate::VerdictStatus;
use std::collections::BTreeMap;

#[test]
fn motifbench_preserves_strict_published_thresholds() {
    let profile = motifbench_1_0();
    let boundary = BTreeMap::from([
        (Box::<str>::from("rmsd"), 2.0),
        (Box::<str>::from("motif_rmsd"), 1.0),
    ]);
    assert_eq!(
        profile.decide_candidate(&boundary).status,
        VerdictStatus::Fail
    );
    assert_eq!(profile.aggregation(), CandidateAggregation::Any);
    assert_eq!(
        profile.metrics()[0].atoms,
        ProfileAtomSet::FullScaffoldCAlpha
    );
    assert_eq!(
        profile.metrics()[1].alignment,
        ProfileAlignment::MeasuredAtoms
    );
}

#[test]
fn ame_keeps_alignment_and_clash_conventions_explicit() {
    let profile = ame_heavy_atom_1_0();
    let passing = BTreeMap::from([
        (Box::<str>::from("catalytic_heavy_atom_rmsd"), 1.499),
        (Box::<str>::from("ligand_backbone_min_distance"), 1.5),
    ]);
    assert_eq!(
        profile.decide_candidate(&passing).status,
        VerdictStatus::Pass
    );
    assert_eq!(
        profile.metrics()[0].alignment,
        ProfileAlignment::CatalyticBackbone
    );
    assert_eq!(
        profile.metrics()[1].alignment,
        ProfileAlignment::PredictionFrame
    );
}

#[test]
fn any_passing_prediction_makes_the_scaffold_pass() {
    let profile = motifbench_1_0();
    let failed = BTreeMap::from([
        (Box::<str>::from("rmsd"), 2.1),
        (Box::<str>::from("motif_rmsd"), 0.5),
    ]);
    let passed = BTreeMap::from([
        (Box::<str>::from("rmsd"), 1.9),
        (Box::<str>::from("motif_rmsd"), 0.9),
    ]);
    let verdict = profile.decide_candidates([&failed, &passed]);
    assert_eq!(verdict.status, VerdictStatus::Pass);
    assert_eq!(verdict.candidates.len(), 2);
}

#[test]
fn an_empty_candidate_set_is_indeterminate() {
    let profile = motifbench_1_0();
    let empty: [&BTreeMap<Box<str>, f64>; 0] = [];
    assert_eq!(
        profile.decide_candidates(empty).status,
        VerdictStatus::Indeterminate
    );
}
