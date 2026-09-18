use super::evaluate_motif;
use crate::{
    AtomSite, Comparison, ComponentRole, ComponentSpec, Constraint, MissingVerdict, Motif,
    NamedConstraint, VerdictProfile, VerdictRule, VerdictStatus,
};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::io::{InputBuffer, ReadOptions};

#[test]
fn orchestration_preserves_each_stage_and_rule() {
    let source = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 O OG SER A 1 0 0 0\nATOM 2 N NE2 HIS A 2 2.8 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("fixture failed: {findings:?}"));
    let motif = Motif::new(
        [
            (
                "ser".into(),
                ComponentSpec::new(ComponentRole::Residue).component("SER"),
            ),
            (
                "his".into(),
                ComponentSpec::new(ComponentRole::Residue).component("HIS"),
            ),
        ],
        [NamedConstraint {
            name: "distance".into(),
            constraint: Constraint::Distance {
                first: AtomSite::new("ser", "OG"),
                second: AtomSite::new("his", "NE2"),
                target: 2.8,
                tolerance: 0.2,
            },
        }],
    )
    .unwrap_or_else(|error| panic!("motif failed: {error}"));
    let profile = VerdictProfile::new(
        "dyad-1.0",
        [VerdictRule {
            metric: "distance".into(),
            comparison: Comparison::Between(2.6, 3.0),
        }],
        MissingVerdict::Indeterminate,
    );
    let report = evaluate_motif(
        &structure,
        &motif,
        None,
        &AnalysisPolicy::default(),
        &profile,
        8,
        crate::MeasurementOptions {
            maximum_alternatives: 8,
            plane_fit: molframe_geom::EigenOptions::standard(),
        },
    )
    .unwrap_or_else(|error| panic!("evaluation failed: {error}"));
    assert!(!report.mapping_ambiguous);
    assert_eq!(report.evaluations.len(), 1);
    let evaluation = &report.evaluations[0];
    assert_eq!(evaluation.verdict.status, VerdictStatus::Pass);
    assert_eq!(evaluation.verdict.outcomes.len(), 1);
    assert_eq!(evaluation.measurements.constraints[0].satisfied, Some(true));
}
