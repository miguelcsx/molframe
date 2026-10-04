use super::*;
use crate::{Code, Diagnostic};

#[test]
fn every_closed_choice_round_trips_through_its_name() {
    for name in Namespace::NAMES {
        assert_eq!(
            name.parse::<Namespace>().expect("parses").to_string(),
            *name
        );
    }
    for name in HydrogenPolicy::NAMES {
        assert_eq!(
            name.parse::<HydrogenPolicy>().expect("parses").name(),
            *name
        );
    }
    for name in RadiiSet::NAMES {
        assert_eq!(name.parse::<RadiiSet>().expect("parses").name(), *name);
    }
    for name in SymmetryPolicy::NAMES {
        assert_eq!(
            name.parse::<SymmetryPolicy>().expect("parses").name(),
            *name
        );
    }
}

#[test]
fn an_underscore_is_read_as_a_hyphen_so_python_spellings_work() {
    assert_eq!(
        "explicit_only".parse::<HydrogenPolicy>(),
        Ok(HydrogenPolicy::ExplicitOnly)
    );
    assert_eq!(
        "amber_united".parse::<RadiiSet>(),
        Ok(RadiiSet::AmberUnited)
    );
    assert_eq!(
        "highest_occupancy_per_atom".parse::<AltlocPolicy>(),
        Ok(AltlocPolicy::HighestOccupancyPerAtom)
    );
    assert_eq!(
        "biological_assembly".parse::<SymmetryPolicy>(),
        Ok(SymmetryPolicy::BiologicalAssembly)
    );
}

#[test]
fn decisions_that_carry_a_value_read_and_print_it() {
    assert_eq!(
        "biological:1".parse::<AssemblyChoice>(),
        Ok(AssemblyChoice::Biological("1".into()))
    );
    assert_eq!(
        "crystal:12.5".parse::<AssemblyChoice>(),
        Ok(AssemblyChoice::Crystal { radius: 12.5 })
    );
    assert_eq!("index:3".parse::<ModelChoice>(), Ok(ModelChoice::Index(3)));
    assert_eq!(
        "label:A".parse::<AltlocPolicy>(),
        Ok(AltlocPolicy::Label("A".into()))
    );
    assert_eq!(
        "explicit:chain A".parse::<AlignmentPolicy>(),
        Ok(AlignmentPolicy::Explicit("chain A".into()))
    );
    assert_eq!(
        "surface:1.4".parse::<ContactDefinition>(),
        Ok(ContactDefinition::SurfaceBased { probe: 1.4 })
    );
    assert_eq!(
        AssemblyChoice::Biological("2".into()).to_string(),
        "biological:2"
    );
    assert_eq!(ModelChoice::Index(7).to_string(), "index:7");
}

#[test]
fn a_word_outside_the_vocabulary_names_the_field_and_the_value() {
    let error = "sometimes"
        .parse::<MissingPolicy>()
        .expect_err("not a policy");
    assert_eq!(error.field, "missing_atoms");
    assert_eq!(error.value, "sometimes");
    for bad in [
        "crystal:0",
        "crystal:-1",
        "crystal:nan",
        "crystal:",
        "biological:",
        "nonsense",
    ] {
        assert!(bad.parse::<AssemblyChoice>().is_err(), "{bad}");
    }
    assert!("index:-1".parse::<ModelChoice>().is_err());
    assert!("distance:-0.1".parse::<ContactDefinition>().is_err());
    assert!("surface:0".parse::<ContactDefinition>().is_err());
}

#[test]
fn a_refused_word_is_a_registered_diagnostic_with_its_field() {
    let diagnostic = Diagnostic::from("x".parse::<Precision>().expect_err("refused"));
    assert_eq!(diagnostic.code(), Code::E6104);
    assert!(diagnostic.code().is_registered());
    assert_eq!(diagnostic.field(), Some("precision"));
}

#[test]
fn polymer_kinds_use_the_shared_vocabulary() {
    use crate::topology::PolymerKind;
    for name in PolymerKind::NAMES {
        assert_eq!(
            name.parse::<PolymerKind>().map(|kind| kind.name()),
            Ok(*name)
        );
    }
    assert_eq!("nucleic_hybrid".parse(), Ok(PolymerKind::NucleicHybrid));
    assert!("lipid".parse::<PolymerKind>().is_err());
}
