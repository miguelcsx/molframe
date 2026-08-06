use super::*;
use crate::contract::policy::{
    AlignmentPolicy, AltlocPolicy, AssemblyChoice, ContactDefinition, EquivalencePolicy,
    Fingerprint, HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy,
    PolicyField, Precision, ProfileId, RadiiSet, SymmetryPolicy,
};

#[test]
fn the_default_policy_is_a_named_profile_rather_than_an_accident() {
    let policy = AnalysisPolicy::default();
    assert_eq!(policy.profile(), Some(ProfileId::DEFAULT));
    assert_eq!(ProfileId::DEFAULT.as_str(), "pdbiox-default-1.0");
}

#[test]
fn the_named_profile_holds_the_values_it_is_documented_to_hold() {
    let policy = AnalysisPolicy::default();
    assert_eq!(policy.assembly, AssemblyChoice::AsymmetricUnit);
    assert_eq!(policy.model, ModelChoice::First);
    assert_eq!(policy.altloc, AltlocPolicy::ConformerConsistent);
    assert_eq!(policy.identifiers, Namespace::Auth);
    assert_eq!(policy.missing_atoms, MissingPolicy::Report);
    assert_eq!(policy.hydrogens, HydrogenPolicy::ExplicitOnly);
    assert_eq!(policy.atom_equivalence, EquivalencePolicy::Ccd);
    assert_eq!(policy.symmetry, SymmetryPolicy::None);
    assert_eq!(policy.alignment, AlignmentPolicy::None);
    assert_eq!(policy.precision, Precision::F64);
    assert_eq!(policy.periodic, PeriodicPolicy::None);
    assert_eq!(policy.vdw_radii, RadiiSet::Bondi);
    assert_eq!(
        policy.contact_def,
        ContactDefinition::DistanceCutoff { tolerance: 0.5 }
    );
}

#[test]
fn per_atom_highest_occupancy_is_never_the_default_and_is_marked_hazardous() {
    assert!(!AnalysisPolicy::default().altloc.is_hazardous());
    assert!(AltlocPolicy::HighestOccupancyPerAtom.is_hazardous());
    assert!(!AltlocPolicy::ConformerConsistent.is_hazardous());
}

#[test]
fn a_policy_altered_anywhere_stops_claiming_to_be_the_profile() {
    for altered in [
        AnalysisPolicy::default().with_identifiers(Namespace::Label),
        AnalysisPolicy::default().with_model(ModelChoice::Index(2)),
        AnalysisPolicy::default().with_altloc(AltlocPolicy::First),
        AnalysisPolicy::default().with_missing_atoms(MissingPolicy::Fail),
        AnalysisPolicy::default().with_assembly(AssemblyChoice::Biological("1".into())),
    ] {
        assert_eq!(altered.profile(), None, "{altered:?}");
    }
}

#[test]
fn a_policy_altered_anywhere_fingerprints_differently() {
    let baseline = AnalysisPolicy::default().fingerprint();
    let variants = [
        AnalysisPolicy::default().with_identifiers(Namespace::Explicit),
        AnalysisPolicy::default().with_model(ModelChoice::All),
        AnalysisPolicy::default().with_altloc(AltlocPolicy::HighestOccupancyPerAtom),
        AnalysisPolicy {
            precision: Precision::F32,
            ..AnalysisPolicy::default()
        },
        AnalysisPolicy {
            contact_def: ContactDefinition::SurfaceBased { probe: 1.4 },
            ..AnalysisPolicy::default()
        },
        AnalysisPolicy {
            vdw_radii: RadiiSet::Charmm,
            ..AnalysisPolicy::default()
        },
    ];
    for variant in variants {
        assert_ne!(variant.fingerprint(), baseline, "{variant:?}");
    }
}

#[test]
fn the_same_policy_fingerprints_the_same_way_every_time() {
    let policy = AnalysisPolicy::default().with_identifiers(Namespace::Explicit);
    assert_eq!(policy.fingerprint(), policy.clone().fingerprint());
    assert_eq!(
        AnalysisPolicy::default().fingerprint(),
        AnalysisPolicy::default().fingerprint()
    );
}

#[test]
fn writing_a_policy_names_every_field_so_the_record_is_complete() {
    let written = AnalysisPolicy::default().to_string();
    for field in PolicyField::ALL {
        assert!(
            written.contains(field.name()),
            "{} is missing",
            field.name()
        );
    }
    assert!(written.contains("pdbiox-default-1.0"));
}

#[test]
fn a_modified_policy_says_so_rather_than_claiming_the_profiles_name() {
    let written = AnalysisPolicy::default()
        .with_identifiers(Namespace::Label)
        .to_string();
    assert!(written.starts_with("(modified from pdbiox-default-1.0)"));
}

#[test]
fn a_fingerprint_of_bytes_depends_on_the_bytes_and_not_on_the_process() {
    assert_eq!(Fingerprint::of(b"abc"), Fingerprint::of(b"abc"));
    assert_ne!(Fingerprint::of(b"abc"), Fingerprint::of(b"abd"));
    assert!(Fingerprint::of(b"abc").to_string().starts_with("fnv1a64:"));
}
