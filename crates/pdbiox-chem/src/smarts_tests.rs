use super::*;
use crate::{ComponentAtom, ComponentBond, ComponentKind, StereoConfiguration};

#[test]
fn carbonyl_and_branch_queries_map_deterministically() {
    let component = acetamide();
    let carbonyl = SmartsPattern::parse("[C;D3](=O)N").expect("valid SMARTS");
    let matches = carbonyl.find_matches(&component);
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].atom_indices.as_ref(), [1, 2, 3]);
}

#[test]
fn aromatic_ring_closure_charge_and_recursive_queries_work() {
    let component = pyridinium();
    assert!(
        SmartsPattern::parse("[n+]1ccccc1")
            .expect("valid")
            .matches(&component)
    );
    assert!(
        SmartsPattern::parse("[r6;x2]")
            .expect("valid")
            .matches(&component)
    );
    assert!(
        SmartsPattern::parse("[$([n+]);R]")
            .expect("valid")
            .matches(&component)
    );
    assert!(
        !SmartsPattern::parse("[n+;H1]")
            .expect("valid")
            .matches(&component)
    );
}

#[test]
fn disconnected_negated_and_stereo_primitives_are_supported() {
    let component = acetamide();
    assert!(
        SmartsPattern::parse("[C;!H1].N")
            .expect("valid")
            .matches(&component)
    );
    assert!(
        SmartsPattern::parse("[C@@]")
            .expect("valid")
            .matches(&component)
    );
    assert!(
        !SmartsPattern::parse("C!=O")
            .expect("valid")
            .matches(&component)
    );
}

#[test]
fn invalid_patterns_report_positions_instead_of_partially_matching() {
    for text in ["", "C(", "C1CC", "[C", "C==O", "[13C]"] {
        assert!(SmartsPattern::parse(text).is_err(), "{text}");
    }
}

fn acetamide() -> Component {
    component(
        &[
            ("CM", Element::CARBON, 0, false, None),
            ("C", Element::CARBON, 0, false, Some(StereoConfiguration::S)),
            ("O", Element::OXYGEN, 0, false, None),
            ("N", Element::NITROGEN, 0, false, None),
        ],
        &[
            ("CM", "C", BondOrder::Single),
            ("C", "O", BondOrder::Double),
            ("C", "N", BondOrder::Single),
        ],
    )
}

fn pyridinium() -> Component {
    component(
        &[
            ("N1", Element::NITROGEN, 1, true, None),
            ("C2", Element::CARBON, 0, true, None),
            ("C3", Element::CARBON, 0, true, None),
            ("C4", Element::CARBON, 0, true, None),
            ("C5", Element::CARBON, 0, true, None),
            ("C6", Element::CARBON, 0, true, None),
        ],
        &[
            ("N1", "C2", BondOrder::Aromatic),
            ("C2", "C3", BondOrder::Aromatic),
            ("C3", "C4", BondOrder::Aromatic),
            ("C4", "C5", BondOrder::Aromatic),
            ("C5", "C6", BondOrder::Aromatic),
            ("C6", "N1", BondOrder::Aromatic),
        ],
    )
}

fn component(
    atoms: &[(&str, Element, i8, bool, Option<StereoConfiguration>)],
    bonds: &[(&str, &str, BondOrder)],
) -> Component {
    Component {
        id: "TST".into(),
        name: "test".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: atoms
            .iter()
            .map(|(name, element, charge, aromatic, stereo)| ComponentAtom {
                name: (*name).into(),
                alternate_name: None,
                element: *element,
                charge: *charge,
                aromatic: *aromatic,
                leaving: false,
                stereo: *stereo,
            })
            .collect::<Vec<_>>()
            .into(),
        bonds: bonds
            .iter()
            .map(|(first, second, order)| ComponentBond {
                atom_a: (*first).into(),
                atom_b: (*second).into(),
                order: *order,
                aromatic: *order == BondOrder::Aromatic,
                stereo: None,
            })
            .collect::<Vec<_>>()
            .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}
