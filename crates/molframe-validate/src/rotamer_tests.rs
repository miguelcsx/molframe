use super::{RotamerDefinition, RotamerOptions, RotamerProfile, rotamer_outliers};
use crate::{ReferenceDistribution, ReferenceLibrary};
use molframe_chem::{Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider};
use molframe_core::contract::DictionaryVersion;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{AnalysisPolicy, BondOrder, Element, Structure};

fn structure() -> Structure {
    let source = "data_r\n\
loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
1 C Q1 CMP A 1 1 0 0\n\
2 C Q2 CMP A 1 0 0 0\n\
3 C Q3 CMP A 1 0 1 0\n\
4 C Q4 CMP A 1 0 1 1\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn provider() -> MemoryProvider {
    let atoms = ["Q1", "Q2", "Q3", "Q4"].map(|name| ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element: Element::CARBON,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    });
    let bonds = [("Q1", "Q2"), ("Q2", "Q3"), ("Q3", "Q4")].map(|(left, right)| ComponentBond {
        atom_a: left.into(),
        atom_b: right.into(),
        order: BondOrder::Single,
        aromatic: false,
        stereo: None,
    });
    MemoryProvider::new(
        DictionaryVersion::new("ccd-test"),
        [Component {
            id: "CMP".into(),
            name: "component".into(),
            kind: ComponentKind::NonPolymer,
            parent: None,
            one_letter_code: None,
            formula: None,
            atoms: atoms.into(),
            bonds: bonds.into(),
            ideal_coordinates: None,
            model_coordinates: None,
        }],
    )
    .expect("component fixture is unique")
}

fn profile() -> RotamerProfile {
    RotamerProfile::new(
        "profile",
        "2026.1",
        [RotamerDefinition {
            component_id: "CMP".into(),
            chi_index: 1,
            atoms: ["Q1".into(), "Q2".into(), "Q3".into(), "Q4".into()],
            distribution: "cmp_chi1".into(),
        }],
    )
    .expect("profile is valid")
}

fn references() -> ReferenceLibrary {
    let distribution =
        ReferenceDistribution::histogram("cmp_chi1", vec![-180.0, 0.0, 180.0], vec![1.0, 1.0])
            .expect("histogram is valid");
    ReferenceLibrary::new("rotamers", "2026.1", [distribution]).expect("set is valid")
}

#[test]
fn explicit_probability_policy_controls_the_outlier_decision() {
    let structure = structure();
    let flagged = rotamer_outliers(
        &structure,
        &provider(),
        &AnalysisPolicy::default(),
        &references(),
        &profile(),
        RotamerOptions {
            minimum_probability: 0.6,
        },
    )
    .expect("validation succeeds");
    assert_eq!(flagged.flags.len(), 1);
    assert_eq!(flagged.flags[0].distribution.as_ref(), "cmp_chi1");
    assert_eq!(flagged.profile_version.as_ref(), "2026.1");
    assert_eq!((flagged.intended, flagged.assessed), (1, 1));

    let accepted = rotamer_outliers(
        &structure,
        &provider(),
        &AnalysisPolicy::default(),
        &references(),
        &profile(),
        RotamerOptions {
            minimum_probability: 0.4,
        },
    )
    .expect("validation succeeds");
    assert!(accepted.flags.is_empty());
}

#[test]
fn non_ccd_connected_profile_path_is_rejected_instead_of_guessed() {
    let invalid = RotamerProfile::new(
        "profile",
        "2026.1",
        [RotamerDefinition {
            component_id: "CMP".into(),
            chi_index: 1,
            atoms: ["Q1".into(), "Q3".into(), "Q2".into(), "Q4".into()],
            distribution: "cmp_chi1".into(),
        }],
    )
    .expect("shape is valid and connectivity is checked against CCD later");
    let result = rotamer_outliers(
        &structure(),
        &provider(),
        &AnalysisPolicy::default(),
        &references(),
        &invalid,
        RotamerOptions {
            minimum_probability: 0.6,
        },
    );
    assert!(matches!(
        result,
        Err(super::RotamerError::InvalidPath { .. })
    ));
}
