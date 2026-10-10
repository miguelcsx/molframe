use super::{MappingOptions, mapped_dockq, mapped_qs_score};
use crate::correspondence::{map_atoms, map_chains, map_residues};
use crate::{DockQOptions, QsOptions};
use molframe_chem::{Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider};
use molframe_core::BondOrder;
use molframe_core::contract::{DictionaryVersion, Namespace, Status};
use molframe_core::element::Element;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_seq::Scoring;
use std::sync::Arc;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(body: &str) -> Structure {
    let input = InputBuffer::from_bytes(format!("{HEADER}{body}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn amino(id: &str, code: u8) -> Component {
    component(id, code, &[], &[])
}

fn component(id: &str, code: u8, atoms: &[&str], bonds: &[(&str, &str)]) -> Component {
    let atoms: Vec<ComponentAtom> = atoms
        .iter()
        .map(|name| ComponentAtom {
            name: (*name).into(),
            alternate_name: None,
            element: Element::from_symbol("C").expect("carbon"),
            charge: 0,
            aromatic: false,
            leaving: false,
            stereo: None,
        })
        .collect();
    let bonds: Vec<ComponentBond> = bonds
        .iter()
        .map(|(first, second)| ComponentBond {
            atom_a: (*first).into(),
            atom_b: (*second).into(),
            order: BondOrder::Single,
            aromatic: false,
            stereo: None,
        })
        .collect();
    Component {
        id: id.into(),
        name: id.into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(code),
        formula: None,
        atoms: Arc::from(atoms),
        bonds: Arc::from(bonds),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn provider() -> MemoryProvider {
    MemoryProvider::new(
        DictionaryVersion::new("test-ccd"),
        [
            amino("ALA", b'A'),
            amino("GLY", b'G'),
            amino("SER", b'S'),
            amino("TRP", b'W'),
            component(
                "FRK",
                b'X',
                &["CB", "CG1", "CG2"],
                &[("CB", "CG1"), ("CB", "CG2")],
            ),
            component(
                "MTH",
                b'M',
                &["CB", "HB1", "HB2", "HB3"],
                &[("CB", "HB1"), ("CB", "HB2"), ("CB", "HB3")],
            ),
        ],
    )
    .expect("component fixtures are unique")
}

fn mapping(provider: &MemoryProvider) -> MappingOptions<'_> {
    MappingOptions {
        provider,
        namespace: Namespace::Label,
        scoring: Scoring::simple(),
        min_identity: 0.5,
        automorphism_limit: 64,
    }
}

// Receptor A (Ala-Gly-Ser) and ligand B (Trp-Ser-Gly) at an interface.
const NATIVE: &str = "\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 2 0 0\n\
ATOM 3 C CA SER A 3 0 2 0\n\
ATOM 4 C CA TRP B 1 1 1 3\n\
ATOM 5 C CA SER B 2 3 1 3\n\
ATOM 6 C CA GLY B 3 1 3 3\n";

// The same complex with chains renamed, listed in the other order, and an
// extra leading serine on the ligand that the native does not have.
const RENAMED: &str = "\
ATOM 1 C CA SER X 1 9 9 9\n\
ATOM 2 C CA TRP X 2 1 1 3\n\
ATOM 3 C CA SER X 3 3 1 3\n\
ATOM 4 C CA GLY X 4 1 3 3\n\
ATOM 5 C CA ALA Y 1 0 0 0\n\
ATOM 6 C CA GLY Y 2 2 0 0\n\
ATOM 7 C CA SER Y 3 0 2 0\n";

fn dockq_options() -> DockQOptions {
    DockQOptions {
        contact_distance: 5.0,
        interface_distance: 10.0,
        ligand_scale: 8.5,
        interface_scale: 1.5,
    }
}

#[test]
fn renamed_reordered_chains_with_an_extra_residue_still_score_one() {
    let native = structure(NATIVE);
    let model = structure(RENAMED);
    let provider = provider();
    let (score, comparison) = mapped_dockq(
        &model,
        &native,
        "A",
        "B",
        &mapping(&provider),
        dockq_options(),
    )
    .expect("the complexes correspond");
    assert!((score.fnat - 1.0).abs() < 1e-12);
    assert!((score.score - 1.0).abs() < 1e-6, "score {}", score.score);
    let pairs: Vec<_> = comparison
        .chains
        .iter()
        .map(|chain| (chain.reference.as_str(), chain.target.as_str()))
        .collect();
    assert_eq!(pairs, [("A", "Y"), ("B", "X")]);
    assert_eq!(comparison.atoms.matches().len(), 6);
}

#[test]
fn a_displaced_ligand_scores_low_through_the_mapping() {
    let native = structure(NATIVE);
    let model = structure(
        "ATOM 1 C CA TRP X 1 51 1 3\n\
ATOM 2 C CA SER X 2 53 1 3\n\
ATOM 3 C CA GLY X 3 51 3 3\n\
ATOM 4 C CA ALA Y 1 0 0 0\n\
ATOM 5 C CA GLY Y 2 2 0 0\n\
ATOM 6 C CA SER Y 3 0 2 0\n",
    );
    let provider = provider();
    let (score, _) = mapped_dockq(
        &model,
        &native,
        "A",
        "B",
        &mapping(&provider),
        dockq_options(),
    )
    .expect("the complexes correspond");
    assert!(score.fnat.abs() < 1e-12);
    assert!(score.score < 0.23, "score {}", score.score);
}

#[test]
fn qs_agrees_across_renamed_chains() {
    let native = structure(NATIVE);
    let model = structure(RENAMED);
    let provider = provider();
    let (score, _) = mapped_qs_score(
        &model,
        &native,
        "A",
        "B",
        &mapping(&provider),
        QsOptions::standard(5.0),
    )
    .expect("the complexes correspond");
    assert!((score - 1.0).abs() < 1e-12);
}

#[test]
fn a_chain_the_model_lacks_is_refused_not_scored() {
    let native = structure(NATIVE);
    let model = structure(
        "ATOM 1 C CA ALA Y 1 0 0 0\n\
ATOM 2 C CA GLY Y 2 2 0 0\n\
ATOM 3 C CA SER Y 3 0 2 0\n",
    );
    let provider = provider();
    let error = mapped_dockq(
        &model,
        &native,
        "A",
        "B",
        &mapping(&provider),
        dockq_options(),
    );
    assert!(error.is_err());
}

#[test]
fn equivalent_atoms_with_swapped_names_are_matched_by_position() {
    let native = structure(
        "ATOM 1 C CB FRK A 1 0 0 0\n\
ATOM 2 C CG1 FRK A 1 1 1 0\n\
ATOM 3 C CG2 FRK A 1 -1 1 0\n",
    );
    let model = structure(
        "ATOM 1 C CB FRK A 1 0 0 0\n\
ATOM 2 C CG1 FRK A 1 -1 1 0\n\
ATOM 3 C CG2 FRK A 1 1 1 0\n",
    );
    let provider = provider();
    let options = mapping(&provider);
    let chains = map_chains(
        &native,
        &model,
        options.provider,
        options.namespace,
        options.scoring,
        options.min_identity,
    )
    .expect("chains map");
    let residues = map_residues(
        &native,
        &model,
        &chains,
        options.provider,
        options.namespace,
        options.scoring,
    )
    .expect("residues map");
    let atoms = map_atoms(&native, &model, &residues, options.provider, 64).expect("atoms map");
    assert_eq!(atoms.swapped_residues, 1);
    let pairs: Vec<_> = atoms
        .mapping
        .matches()
        .iter()
        .map(|pair| (pair.reference, pair.model))
        .collect();
    assert_eq!(pairs, [(0, 0), (1, 2), (2, 1)]);
}

#[test]
fn equivalences_are_counted_over_the_atoms_both_residues_have() {
    // The dictionary's three hydrogens give six equivalent mappings; neither
    // structure has any of them, so the bound of two is never approached.
    let body = "ATOM 1 C CB MTH A 1 0 0 0\n";
    let native = structure(body);
    let model = structure(body);
    let provider = provider();
    let options = mapping(&provider);
    let chains = map_chains(
        &native,
        &model,
        options.provider,
        options.namespace,
        options.scoring,
        options.min_identity,
    )
    .expect("chains map");
    let residues = map_residues(
        &native,
        &model,
        &chains,
        options.provider,
        options.namespace,
        options.scoring,
    )
    .expect("residues map");
    let atoms = map_atoms(&native, &model, &residues, options.provider, 2).expect("atoms map");
    assert_eq!(atoms.status, Status::Complete);
    assert_eq!(atoms.mapping.matches().len(), 1);
}

#[test]
fn equal_cost_ccd_equivalent_mappings_are_reported_as_ambiguous() {
    let native = structure(
        "ATOM 1 C CB FRK A 1 0 0 0\n\
ATOM 2 C CG1 FRK A 1 1 1 0\n\
ATOM 3 C CG2 FRK A 1 -1 1 0\n",
    );
    let model = structure(
        "ATOM 1 C CB FRK A 1 0 0 0\n\
ATOM 2 C CG1 FRK A 1 0 1 0\n\
ATOM 3 C CG2 FRK A 1 0 1 0\n",
    );
    let provider = provider();
    let options = mapping(&provider);
    let chains = map_chains(
        &native,
        &model,
        options.provider,
        options.namespace,
        options.scoring,
        options.min_identity,
    )
    .expect("chains map");
    let residues = map_residues(
        &native,
        &model,
        &chains,
        options.provider,
        options.namespace,
        options.scoring,
    )
    .expect("residues map");
    let atoms = map_atoms(&native, &model, &residues, options.provider, 64).expect("atoms map");

    assert_eq!(atoms.status, Status::Ambiguous);
    assert_eq!(atoms.alternatives.len(), 1);
    assert_eq!(atoms.alternatives[0].mappings.len(), 2);
    assert_eq!(atoms.mapping.matches().len(), 3);
}

// A homodimer: both chains read ALA-GLY-SER, so sequence cannot say which model
// chain is which native chain.
const HOMODIMER_NATIVE: &str = "\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 2 0 0\n\
ATOM 3 C CA SER A 3 0 2 0\n\
ATOM 4 C CA ALA B 1 1 1 3\n\
ATOM 5 C CA GLY B 2 3 1 3\n\
ATOM 6 C CA SER B 3 1 3 3\n";

// The same dimer with the two chains exchanged: X sits where B was, Y where A
// was, and X is listed first.
const HOMODIMER_SWAPPED: &str = "\
ATOM 1 C CA ALA X 1 1 1 3\n\
ATOM 2 C CA GLY X 2 3 1 3\n\
ATOM 3 C CA SER X 3 1 3 3\n\
ATOM 4 C CA ALA Y 1 0 0 0\n\
ATOM 5 C CA GLY Y 2 2 0 0\n\
ATOM 6 C CA SER Y 3 0 2 0\n";

#[test]
fn sequence_alone_would_pair_the_swapped_homodimer_chains_in_label_order() {
    let native = structure(HOMODIMER_NATIVE);
    let model = structure(HOMODIMER_SWAPPED);
    let provider = provider();
    let options = mapping(&provider);
    let chains = map_chains(
        &native,
        &model,
        options.provider,
        options.namespace,
        options.scoring,
        options.min_identity,
    )
    .expect("sequences map");
    let pairs: Vec<_> = chains
        .iter()
        .map(|chain| (chain.reference.as_str(), chain.target.as_str()))
        .collect();
    assert_eq!(pairs, [("A", "X"), ("B", "Y")]);
}

#[test]
fn a_swapped_homodimer_scores_one_under_the_assignment_that_fits() {
    let native = structure(HOMODIMER_NATIVE);
    let model = structure(HOMODIMER_SWAPPED);
    let provider = provider();
    let (score, comparison) = mapped_dockq(
        &model,
        &native,
        "A",
        "B",
        &mapping(&provider),
        dockq_options(),
    )
    .expect("the complexes correspond");
    assert!((score.score - 1.0).abs() < 1e-6, "score {}", score.score);
    assert_eq!(comparison.chain_assignments_tried, 2);
    let pairs: Vec<_> = comparison
        .chains
        .iter()
        .map(|chain| (chain.reference.as_str(), chain.target.as_str()))
        .collect();
    assert_eq!(pairs, [("A", "Y"), ("B", "X")]);
}

#[test]
fn qs_reports_how_many_tied_assignments_it_compared() {
    let native = structure(HOMODIMER_NATIVE);
    let model = structure(HOMODIMER_SWAPPED);
    let provider = provider();
    let (score, comparison) = mapped_qs_score(
        &model,
        &native,
        "A",
        "B",
        &mapping(&provider),
        QsOptions::standard(5.0),
    )
    .expect("the complexes correspond");
    assert!((score - 1.0).abs() < 1e-12);
    assert_eq!(comparison.chain_assignments_tried, 2);
}

#[test]
fn distinguishable_chains_are_assigned_once() {
    let native = structure(NATIVE);
    let model = structure(RENAMED);
    let provider = provider();
    let (_, comparison) = mapped_dockq(
        &model,
        &native,
        "A",
        "B",
        &mapping(&provider),
        dockq_options(),
    )
    .expect("the complexes correspond");
    assert_eq!(comparison.chain_assignments_tried, 1);
}
