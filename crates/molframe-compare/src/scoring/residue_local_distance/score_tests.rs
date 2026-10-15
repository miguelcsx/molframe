use super::{ResidueCorrespondence, ResidueLddtOptions, SymmetryPolicy, residue_lddt};
use crate::pair_residues_by_identifier;
use molframe_chem::{Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider};
use molframe_core::contract::{DictionaryVersion, Namespace};
use molframe_core::element::Element;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_core::{BondOrder, ExecutionContext};
use std::fmt::Write;

type Atom = (&'static str, f32, f32, f32);

const ASP: [Atom; 8] = [
    ("N", 0.0, 0.0, 0.0),
    ("CA", 1.4, 0.4, 0.1),
    ("C", 2.5, -0.5, 0.9),
    ("O", 2.4, -1.7, 1.1),
    ("CB", 1.6, 1.9, 0.4),
    ("CG", 2.9, 2.5, -0.2),
    ("OD1", 3.1, 3.7, -0.4),
    ("OD2", 3.9, 1.8, -0.6),
];
const ALA: [Atom; 5] = [
    ("N", 0.0, 0.0, 0.0),
    ("CA", 1.4, 0.4, 0.1),
    ("C", 2.5, -0.5, 0.9),
    ("O", 2.4, -1.7, 1.1),
    ("CB", 1.6, 1.9, 0.4),
];
const VAL: [Atom; 7] = [
    ("N", 0.0, 0.0, 0.0),
    ("CA", 1.4, 0.4, 0.1),
    ("C", 2.5, -0.5, 0.9),
    ("O", 2.4, -1.7, 1.1),
    ("CB", 1.6, 1.9, 0.4),
    ("CG1", 0.9, 3.1, 0.2),
    ("CG2", 2.9, 2.3, 1.5),
];

fn element(name: &str) -> Element {
    let symbol = &name[..1];
    match Element::from_symbol(symbol) {
        Some(element) => element,
        None => panic!("unknown element {symbol}"),
    }
}

fn component(id: &str, atoms: &[Atom], bonds: &[(&str, &str, BondOrder)]) -> Component {
    Component {
        id: id.into(),
        name: id.into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(b'X'),
        formula: None,
        atoms: atoms
            .iter()
            .map(|(name, ..)| ComponentAtom {
                name: (*name).into(),
                alternate_name: None,
                element: element(name),
                charge: 0,
                aromatic: false,
                leaving: false,
                stereo: None,
            })
            .collect(),
        bonds: bonds
            .iter()
            .map(|(first, second, order)| ComponentBond {
                atom_a: (*first).into(),
                atom_b: (*second).into(),
                order: *order,
                aromatic: false,
                stereo: None,
            })
            .collect(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn provider() -> MemoryProvider {
    use BondOrder::{Double, Single};
    let backbone = [
        ("N", "CA", Single),
        ("CA", "C", Single),
        ("C", "O", Double),
        ("CA", "CB", Single),
    ];
    let with = |extra: &[(&'static str, &'static str, BondOrder)]| {
        let mut bonds = backbone.to_vec();
        bonds.extend_from_slice(extra);
        bonds
    };
    MemoryProvider::new(
        DictionaryVersion::new("test-ccd"),
        [
            component(
                "ASP",
                &ASP,
                &with(&[
                    ("CB", "CG", Single),
                    ("CG", "OD1", Double),
                    ("CG", "OD2", Single),
                ]),
            ),
            component("ALA", &ALA, &with(&[])),
            component(
                "VAL",
                &VAL,
                &with(&[("CB", "CG1", Single), ("CB", "CG2", Single)]),
            ),
        ],
    )
    .expect("component fixtures are unique")
}

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

/// Residues placed 3.8 Å apart along x; `edit` may drop or move atoms.
fn build(
    residues: &[(&str, &[Atom])],
    edit: impl Fn(usize, &str, [f32; 3]) -> Option<[f32; 3]>,
) -> Structure {
    let mut body = String::new();
    let mut serial = 0;
    for (index, (name, atoms)) in residues.iter().enumerate() {
        for (atom, x, y, z) in *atoms {
            let moved = [x + 3.8 * shift(index), *y, *z];
            let Some(position) = edit(index, atom, moved) else {
                continue;
            };
            serial += 1;
            let _ = writeln!(
                body,
                "ATOM {serial} {} {atom} {name} A {} {} {} {}",
                &atom[..1],
                index + 1,
                position[0],
                position[1],
                position[2],
            );
        }
    }
    let input = InputBuffer::from_bytes(format!("{HEADER}{body}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

const CHAIN: [(&str, &[Atom]); 3] = [("ASP", &ASP), ("ALA", &ALA), ("VAL", &VAL)];

fn intact(residues: &[(&str, &[Atom])]) -> Structure {
    build(residues, |_, _, position| Some(position))
}

fn shift(index: usize) -> f32 {
    let Ok(small) = u8::try_from(index) else {
        panic!("fixture index out of range");
    };
    f32::from(small)
}

fn score(
    reference: &Structure,
    model: &Structure,
    options: &ResidueLddtOptions,
) -> super::ResidueLddt {
    let provider = provider();
    let pairs = match pair_residues_by_identifier(reference, model, &provider, Namespace::Label) {
        Ok(pairs) => pairs,
        Err(error) => panic!("pairing failed: {error:?}"),
    };
    match residue_lddt(
        reference,
        model,
        &pairs,
        &provider,
        options,
        &ExecutionContext::default(),
    ) {
        Ok(result) => result,
        Err(error) => panic!("score failed: {error}"),
    }
}

fn without_symmetry() -> ResidueLddtOptions {
    ResidueLddtOptions {
        symmetry: SymmetryPolicy::Disabled,
        ..ResidueLddtOptions::standard()
    }
}

#[test]
fn a_model_identical_to_the_reference_scores_one() {
    let reference = intact(&CHAIN);
    let result = score(&reference, &reference, &ResidueLddtOptions::standard());
    assert!((result.score - 1.0).abs() < 1e-12);
    assert!(result.considered > 0);
    assert_eq!(result.renamed_residues, 0);
    assert!(result.residues.iter().all(|r| r.score == Some(1.0)));
}

fn swap_pair(first: &'static str, second: &'static str, residue: usize) -> Structure {
    let find = |name: &str| {
        let atoms: &[Atom] = CHAIN[residue].1;
        let (_, x, y, z) = atoms.iter().find(|atom| atom.0 == name).expect("atom");
        [x + 3.8 * shift(residue), *y, *z]
    };
    let (a, b) = (find(first), find(second));
    build(&CHAIN, move |index, name, position| {
        if index != residue {
            Some(position)
        } else if name == first {
            Some(b)
        } else if name == second {
            Some(a)
        } else {
            Some(position)
        }
    })
}

#[test]
fn swapped_carboxylate_oxygens_are_not_penalised() {
    let reference = intact(&CHAIN);
    let model = swap_pair("OD1", "OD2", 0);
    let resolved = score(&reference, &model, &ResidueLddtOptions::standard());
    assert!((resolved.score - 1.0).abs() < 1e-12, "{}", resolved.score);
    assert_eq!(resolved.renamed_residues, 1);
    let literal = score(&reference, &model, &without_symmetry());
    assert!(literal.score < 1.0, "{}", literal.score);
}

#[test]
fn the_two_methyls_of_valine_are_not_interchangeable() {
    let reference = intact(&CHAIN);
    let model = swap_pair("CG1", "CG2", 2);
    let result = score(&reference, &model, &ResidueLddtOptions::standard());
    assert!(result.score < 1.0, "{}", result.score);
    assert_eq!(result.renamed_residues, 0);
}

#[test]
fn a_missing_model_atom_lowers_the_score_of_its_residue() {
    let reference = intact(&CHAIN);
    let model = build(&CHAIN, |index, name, position| {
        (index != 1 || name != "CB").then_some(position)
    });
    let result = score(&reference, &model, &ResidueLddtOptions::standard());
    assert!(result.score < 1.0);
    assert!(result.residues[1].score.is_some_and(|value| value < 1.0));
    assert!(result.residues[2].score.is_some_and(|value| value < 1.0));
}

#[test]
fn a_missing_residue_counts_against_the_model() {
    let reference = intact(&CHAIN);
    let model = intact(&CHAIN[..2]);
    let result = score(&reference, &model, &ResidueLddtOptions::standard());
    assert_eq!(result.residues[2].model, None);
    assert_eq!(result.residues[2].score, Some(0.0));
    assert!(result.score < 1.0);
}

#[test]
fn distances_inside_one_residue_are_never_tested() {
    let reference = intact(&CHAIN[..1]);
    let model = build(&CHAIN[..1], |_, name, position| {
        (name == "OD1")
            .then_some([position[0] + 9.0, position[1], position[2]])
            .or(Some(position))
    });
    let result = score(&reference, &model, &ResidueLddtOptions::standard());
    assert_eq!(result.considered, 0);
    let strict = ResidueLddtOptions {
        empty_policy: crate::EmptyLddtPolicy::Error,
        ..ResidueLddtOptions::standard()
    };
    let provider = provider();
    let pairs = [ResidueCorrespondence::from((
        reference.data().residues().next().expect("residue").index(),
        model.data().residues().next().expect("residue").index(),
    ))];
    assert!(
        residue_lddt(
            &reference,
            &model,
            &pairs,
            &provider,
            &strict,
            &ExecutionContext::default()
        )
        .is_err()
    );
}

#[test]
fn sequence_separation_removes_pairs_between_nearby_residues() {
    let reference = intact(&CHAIN);
    let near = score(&reference, &reference, &ResidueLddtOptions::standard());
    let far = score(
        &reference,
        &reference,
        &ResidueLddtOptions {
            sequence_separation: 1,
            ..ResidueLddtOptions::standard()
        },
    );
    assert!(far.considered > 0 && far.considered < near.considered);
}

#[test]
fn identifier_pairing_leaves_an_absent_residue_unpaired() {
    let reference = intact(&CHAIN);
    let model = intact(&CHAIN[..2]);
    let pairs = pair_residues_by_identifier(&reference, &model, &provider(), Namespace::Label)
        .expect("pairing");
    assert_eq!(pairs.len(), 3);
    assert!(pairs[0].model.is_some() && pairs[2].model.is_none());
}

#[test]
fn a_coordinate_scramble_scores_far_below_a_faithful_model() {
    let reference = intact(&CHAIN);
    let model = build(&CHAIN, |index, _, position| {
        Some([
            position[0] * 2.0 + shift(index),
            -position[1],
            position[2] + 3.0,
        ])
    });
    let result = score(&reference, &model, &ResidueLddtOptions::standard());
    assert!(result.score < 0.5, "{}", result.score);
}

#[test]
fn the_governed_score_reports_coverage_and_every_convention() {
    let reference = intact(&CHAIN);
    let model = intact(&CHAIN[..2]);
    let provider = provider();
    let pairs = pair_residues_by_identifier(&reference, &model, &provider, Namespace::Label)
        .expect("pairing");
    let policy = molframe_core::contract::AnalysisPolicy::default();
    let analysis = crate::governed_residue_lddt(
        &reference,
        &model,
        &pairs,
        &provider,
        &ResidueLddtOptions::standard(),
        &policy,
        &ExecutionContext::default(),
    );
    let Ok(analysis) = analysis else {
        panic!("governed score failed");
    };
    assert_eq!(analysis.coverage.intended, 3);
    assert_eq!(analysis.coverage.used, 2);
    let parameters = &analysis.provenance.parameters;
    for key in [
        "symmetry",
        "missing_model_atoms",
        "pair_counting",
        "tolerances",
    ] {
        assert!(parameters.contains_key(key), "{key}");
    }
}
