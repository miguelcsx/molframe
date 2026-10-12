use super::*;
use crate::bond::{BondProvenance, BondRecord, BondTableBuilder};
use crate::structure::fixture;

fn options(tolerance: f32) -> StructureDifferenceOptions {
    StructureDifferenceOptions {
        coordinate_tolerance: tolerance,
    }
}

/// A structure with one bond whose connection class is interned locally.
///
/// `filler` is interned first, when given, so the connection class lands in a
/// known slot. Each call clones a fixed fixture dictionary and interns from
/// there, so the identifier a class receives depends on what came before it —
/// which is exactly the property two independently built structures have.
///
/// The class names must not be canonical strings: those have fixed identifiers
/// shared by every dictionary, which would make these tests pass by accident.
fn classified(filler: Option<&str>, class: &str) -> Structure {
    let mut data = fixture::sample().data().clone();
    if let Some(filler) = filler {
        let _ = data.dictionary.intern(filler);
    }
    let Ok(kind) = data.dictionary.intern(class) else {
        panic!("interning `{class}` failed")
    };
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: crate::AtomIndex::new(0),
        atom_b: crate::AtomIndex::new(1),
        order: crate::BondOrder::Single,
        provenance: BondProvenance::File,
        kind: crate::optional::OptionalSymbol::some(kind),
    });
    data.bonds = bonds.finish();
    Structure::new(data)
}

fn class_id(structure: &Structure) -> Option<crate::optional::OptionalSymbol> {
    structure
        .data()
        .bonds
        .get(crate::BondIndex::new(0))
        .map(|bond| bond.kind)
}

#[test]
fn identical_snapshots_have_no_semantic_difference() {
    let structure = fixture::sample();
    let result = structure_difference(&structure, &structure, options(0.0));
    assert!(matches!(result, Ok(difference) if difference.is_empty()));
}

/// The same connection class, numbered differently by two dictionaries, is not a
/// difference.
///
/// A connection class is an interned identifier, so where a class lands depends
/// on what its dictionary saw first. Comparing the raw values would report "these
/// structures differ" for two files that record exactly the same chemistry.
#[test]
fn the_same_connection_class_in_two_dictionaries_is_not_a_difference() {
    // The same string, but the right-hand dictionary interned a filler first, so
    // it lands in a later dynamic slot.
    let left = classified(None, "source-class-alpha");
    let right = classified(Some("source-class-filler"), "source-class-alpha");
    assert_ne!(
        class_id(&left),
        class_id(&right),
        "the fixture must number the same class differently"
    );

    let Ok(difference) = structure_difference(&left, &right, options(0.0)) else {
        panic!("a semantic comparison of two valid structures must succeed");
    };
    assert_eq!(
        difference.changed_bonds, 0,
        "the same connection class in two dictionaries is not a change"
    );
    assert!(difference.is_empty(), "{difference:?}");
}

/// Two different connection classes that happen to share a numeric identifier are
/// a difference.
///
/// Each string is the first non-canonical entry in its own dictionary, so both
/// receive the same slot. A raw comparison would call them equal and report two
/// files with genuinely different connectivity as identical, which is the more
/// dangerous direction: it hides a change.
#[test]
fn different_connection_classes_sharing_one_identifier_are_a_difference() {
    let left = classified(None, "source-class-alpha");
    let right = classified(None, "source-class-beta");
    assert_eq!(
        class_id(&left),
        class_id(&right),
        "the fixture must collide the two classes into one identifier"
    );

    let Ok(difference) = structure_difference(&left, &right, options(0.0)) else {
        panic!("a semantic comparison of two valid structures must succeed");
    };
    assert_eq!(
        difference.changed_bonds, 1,
        "two different connection classes are different connectivity"
    );
}

#[test]
fn metadata_and_coordinates_are_reported_independently() {
    let left = fixture::sample();
    let mut data = left.data().clone();
    data.entry.title = Some("changed".into());
    let right = Structure::new(data);
    let result = structure_difference(&left, &right, options(0.1));
    let Ok(difference) = result else {
        panic!("valid comparison failed")
    };
    assert!(!difference.is_empty(), "{difference:?}");
}

#[test]
fn coordinate_tolerance_must_be_explicitly_valid() {
    let structure = fixture::sample();
    assert_eq!(
        structure_difference(&structure, &structure, options(f32::NAN)),
        Err(DifferenceError::InvalidCoordinateTolerance)
    );
    assert_eq!(
        structure_difference(&structure, &structure, options(-0.1)),
        Err(DifferenceError::InvalidCoordinateTolerance)
    );
}
