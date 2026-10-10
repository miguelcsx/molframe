use super::super::*;
use crate::{ComponentAtom, ComponentBond, MemoryProvider, StereoConfiguration};
use molframe_core::contract::DictionaryVersion;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{BondOrder, Element};

const HEADER: &str = "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn comp(id: &str, kind: ComponentKind, atoms: &[(&str, Element)]) -> Component {
    Component {
        id: id.into(),
        name: id.into(),
        kind,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: atoms
            .iter()
            .map(|(name, element)| ComponentAtom {
                name: (*name).into(),
                alternate_name: None,
                element: *element,
                charge: 0,
                aromatic: false,
                leaving: false,
                stereo: None,
            })
            .collect::<Vec<_>>()
            .into(),
        bonds: Vec::<ComponentBond>::new().into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn gly() -> Component {
    comp(
        "GLY",
        ComponentKind::AminoAcid,
        &[
            ("C", Element::CARBON),
            ("N", Element::NITROGEN),
            ("N2", Element::NITROGEN),
        ],
    )
}

fn run(body: &str, components: Vec<Component>, rules: Vec<PolymerLinkRule>) -> ChemistryReport {
    let text = format!("{HEADER}{body}");
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider = MemoryProvider::new(DictionaryVersion::new("t"), components).expect("unique");
    apply_component_chemistry(
        &structure,
        &provider,
        PolymerLinkPolicy::explicit(2.0, rules),
    )
    .expect("annotation")
}

fn polymeric(report: &ChemistryReport) -> usize {
    report
        .structure
        .data()
        .bonds
        .iter()
        .filter(|bond| bond.order == BondOrder::Polymeric)
        .count()
}

fn peptide_rule() -> PolymerLinkRule {
    PolymerLinkRule::new(ComponentKind::AminoAcid, "C", ComponentKind::AminoAcid, "N")
}

fn alt_rule() -> PolymerLinkRule {
    PolymerLinkRule::new(
        ComponentKind::AminoAcid,
        "C",
        ComponentKind::AminoAcid,
        "N2",
    )
}

fn reason(report: &ChemistryReport) -> Option<String> {
    report
        .findings
        .iter()
        .find(|finding| finding.code() == Code::W3302)
        .and_then(|finding| {
            finding
                .context()
                .iter()
                .find(|item| item.label() == "reason")
                .map(|item| item.value().to_string())
        })
}

#[test]
fn a_later_matching_rule_is_evaluated_when_an_earlier_rule_has_no_atoms() {
    let body = "ATOM 1 C C . GLY A 1 0 0 0\nATOM 2 N N . GLY A 2 1.4 0 0\n";
    let missing =
        PolymerLinkRule::new(ComponentKind::AminoAcid, "X", ComponentKind::AminoAcid, "Y");
    let report = run(body, vec![gly()], vec![missing, peptide_rule()]);
    assert_eq!(polymeric(&report), 1);
}

#[test]
fn altloc_conformers_link_only_to_their_own_conformer() {
    let body = "ATOM 1 C C A GLY A 1 0 0 0\nATOM 2 C C B GLY A 1 0 0.3 0\n\
ATOM 3 N N A GLY A 2 1.4 0 0\nATOM 4 N N B GLY A 2 1.4 0.3 0\n";
    let report = run(body, vec![gly()], vec![peptide_rule()]);
    // A-A and B-B only; the cross pairs are 0.3 A apart in y and compatible
    // by distance alone but belong to different conformers.
    assert_eq!(polymeric(&report), 2);
    let pairs: Vec<(usize, usize)> = report
        .structure
        .data()
        .bonds
        .iter()
        .map(|bond| (bond.atom_a.as_usize(), bond.atom_b.as_usize()))
        .collect();
    assert!(
        pairs.contains(&(0, 2)) && pairs.contains(&(1, 3)),
        "{pairs:?}"
    );
}

#[test]
fn a_too_short_contact_is_not_a_bond_and_is_diagnosed() {
    let body = "ATOM 1 C C . GLY A 1 0 0 0\nATOM 2 N N . GLY A 2 0.5 0 0\n";
    let report = run(body, vec![gly()], vec![peptide_rule()]);
    assert_eq!(polymeric(&report), 0);
    assert_eq!(reason(&report).as_deref(), Some("too short"));
}

#[test]
fn a_distant_contact_is_diagnosed_as_a_gap() {
    let body = "ATOM 1 C C . GLY A 1 0 0 0\nATOM 2 N N . GLY A 2 1.9 0 0\n";
    let report = run(body, vec![gly()], vec![peptide_rule()]);
    // 1.9 A exceeds the C-N window even though the policy allows 2.0 A.
    assert_eq!(polymeric(&report), 0);
    assert_eq!(reason(&report).as_deref(), Some("too far"));
}

#[test]
fn equidistant_competing_candidates_are_reported_not_guessed() {
    let body = "ATOM 1 C C . GLY A 1 0 0 0\nATOM 2 N N . GLY A 2 1.4 0 0\n\
ATOM 3 N N2 . GLY A 2 -1.4 0 0\n";
    let report = run(body, vec![gly()], vec![peptide_rule(), alt_rule()]);
    assert_eq!(polymeric(&report), 0);
    assert_eq!(
        reason(&report).as_deref(),
        Some("ambiguous attachment candidates")
    );
}

#[test]
fn the_shortest_candidate_wins_an_attachment_atom() {
    let body = "ATOM 1 C C . GLY A 1 0 0 0\nATOM 2 N N . GLY A 2 1.3 0 0\n\
ATOM 3 N N2 . GLY A 2 -1.5 0 0\n";
    let report = run(body, vec![gly()], vec![peptide_rule(), alt_rule()]);
    assert_eq!(polymeric(&report), 1);
    assert!(reason(&report).is_none());
}

#[test]
fn a_protein_chain_with_water_and_ligands_is_still_protein() {
    let hoh = comp("HOH", ComponentKind::Solvent, &[("O", Element::OXYGEN)]);
    let lig = comp("LIG", ComponentKind::NonPolymer, &[("O", Element::OXYGEN)]);
    let body = "ATOM 1 C C . GLY A 1 0 0 0\nATOM 2 C C . GLY A 2 5 0 0\n\
ATOM 3 C C . GLY A 3 10 0 0\nHETATM 4 O O . HOH A 4 15 0 0\n\
HETATM 5 O O . LIG A 5 20 0 0\n";
    let report = run(body, vec![gly(), hoh, lig], vec![peptide_rule()]);
    assert_eq!(
        report
            .structure
            .chain(molframe_core::ChainIndex::new(0))
            .map(molframe_core::structure::ChainRef::polymer_kind),
        Some(PolymerKind::Protein)
    );
}

#[test]
fn formal_charges_are_recorded_as_component_defaults() {
    let report = run(
        "ATOM 1 C C . GLY A 1 0 0 0\n",
        vec![gly()],
        vec![peptide_rule()],
    );
    assert_eq!(
        report.provenance.charge_source,
        FormalChargeSource::ComponentDefault
    );
}

#[test]
fn stereo_uses_only_atoms_of_the_centres_own_conformer() {
    let ligand = Component {
        atoms: vec![
            ("CTR", Element::CARBON, Some(StereoConfiguration::R)),
            ("A", Element::NITROGEN, None),
            ("B", Element::OXYGEN, None),
            ("C", Element::SULFUR, None),
        ]
        .into_iter()
        .map(|(name, element, stereo)| ComponentAtom {
            name: name.into(),
            alternate_name: None,
            element,
            charge: 0,
            aromatic: false,
            leaving: false,
            stereo,
        })
        .collect::<Vec<_>>()
        .into(),
        bonds: ["A", "B", "C"]
            .into_iter()
            .map(|other| ComponentBond {
                atom_a: "CTR".into(),
                atom_b: other.into(),
                order: BondOrder::Single,
                aromatic: false,
                stereo: None,
            })
            .collect::<Vec<_>>()
            .into(),
        ideal_coordinates: Some(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
            .into(),
        ),
        ..comp("LIG", ComponentKind::NonPolymer, &[])
    };
    // Conformer A places N at +x, conformer B at -x. The centre is conformer B.
    // Reading the first "A" atom (conformer A) would give S; conformer B is R.
    let body = "HETATM 1 C CTR B LIG A 1 0 0 0\nHETATM 2 N A A LIG A 1 1 0 0\n\
HETATM 3 N A B LIG A 1 -1 0 0\nHETATM 4 O B . LIG A 1 0 1 0\n\
HETATM 5 S C . LIG A 1 0 0 -1\n";
    let report = run(body, vec![ligand], vec![peptide_rule()]);
    let Some(molframe_core::AtomAnnotation::Symbol(stereo)) = report
        .structure
        .annotations()
        .get(molframe_core::STEREO_CONFIGURATION_ANNOTATION)
    else {
        panic!("stereo annotation absent")
    };
    let Some((symbol, molframe_core::Presence::Present)) = stereo.get(0) else {
        panic!("centre configuration absent")
    };
    assert_eq!(report.structure.resolve(symbol), Some("R"));
}
