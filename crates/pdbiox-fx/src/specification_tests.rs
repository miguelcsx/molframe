use super::{AtomSite, ComponentRole, ComponentSpec, Constraint, Motif, NamedConstraint};

#[test]
fn catalytic_triad_is_a_valid_immutable_specification() {
    let motif = Motif::new(
        [
            (
                "ser".into(),
                ComponentSpec::new(ComponentRole::Residue)
                    .component("SER")
                    .require("OG"),
            ),
            (
                "his".into(),
                ComponentSpec::new(ComponentRole::Residue)
                    .component("HIS")
                    .require("NE2"),
            ),
            (
                "asp".into(),
                ComponentSpec::new(ComponentRole::Residue)
                    .component("ASP")
                    .equivalent(["OD1", "OD2"]),
            ),
        ],
        [NamedConstraint {
            name: "ser-his".into(),
            constraint: Constraint::Distance {
                first: AtomSite::new("ser", "OG"),
                second: AtomSite::new("his", "NE2"),
                target: 2.8,
                tolerance: 0.4,
            },
        }],
    )
    .unwrap_or_else(|error| panic!("motif failed: {error}"));
    assert_eq!(motif.components().len(), 3);
    assert_eq!(motif.constraints().len(), 1);
}

#[test]
fn a_constraint_cannot_reference_an_undeclared_component() {
    let result = Motif::new(
        [("ser".into(), ComponentSpec::new(ComponentRole::Residue))],
        [NamedConstraint {
            name: "bad".into(),
            constraint: Constraint::StericExclusion {
                first: AtomSite::new("ser", "OG"),
                second: AtomSite::new("missing", "C"),
                min_distance: 2.0,
            },
        }],
    );
    assert!(result.is_err());
}
