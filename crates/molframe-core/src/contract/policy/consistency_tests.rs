use super::*;

fn policy(assembly: AssemblyChoice, symmetry: SymmetryPolicy) -> AnalysisPolicy {
    AnalysisPolicy {
        assembly,
        symmetry,
        ..AnalysisPolicy::default()
    }
}

#[test]
fn the_defaults_and_each_matching_pair_are_consistent() {
    assert!(AnalysisPolicy::default().check_consistency().is_ok());
    let biological = AssemblyChoice::Biological("1".into());
    let crystal = AssemblyChoice::Crystal { radius: 8.0 };
    for (assembly, symmetry) in [
        (biological.clone(), SymmetryPolicy::None),
        (biological, SymmetryPolicy::BiologicalAssembly),
        (crystal.clone(), SymmetryPolicy::None),
        (crystal, SymmetryPolicy::Crystallographic),
    ] {
        assert!(policy(assembly, symmetry).check_consistency().is_ok());
    }
}

#[test]
fn a_pair_that_describes_no_single_system_names_both_fields() {
    for (assembly, symmetry) in [
        (
            AssemblyChoice::AsymmetricUnit,
            SymmetryPolicy::Crystallographic,
        ),
        (
            AssemblyChoice::AsymmetricUnit,
            SymmetryPolicy::BiologicalAssembly,
        ),
        (
            AssemblyChoice::Biological("1".into()),
            SymmetryPolicy::Crystallographic,
        ),
        (
            AssemblyChoice::Crystal { radius: 8.0 },
            SymmetryPolicy::BiologicalAssembly,
        ),
    ] {
        let Err(finding) = policy(assembly, symmetry).check_consistency() else {
            panic!("a contradictory pair must be refused");
        };
        assert_eq!(finding.code(), Code::E6004);
        let context = finding.context();
        assert!(
            context
                .iter()
                .any(|item| item.label() == "fields" && item.value() == "assembly, symmetry")
        );
    }
}

#[test]
fn a_crystal_radius_must_be_positive_and_finite() {
    for radius in [0.0, -1.0, f32::NAN, f32::INFINITY] {
        let Err(finding) =
            policy(AssemblyChoice::Crystal { radius }, SymmetryPolicy::None).check_consistency()
        else {
            panic!("radius {radius} must be refused");
        };
        assert_eq!(finding.code(), Code::E6016);
    }
}
