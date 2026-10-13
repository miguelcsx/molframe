use super::*;
use crate::Component;
use crate::{ComponentAtom, ComponentBond, LazyCifProvider, MemoryProvider, StereoConfiguration};
use molframe_core::bond::{BondProvenance, BondRecord};
use molframe_core::contract::DictionaryVersion;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{BondOrder, Element};
use std::fmt::Write as _;

const STRUCTURE: &str = "data_gly\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N GLY A 1 0.0 0 0\nATOM 2 C CA GLY A 1 1.4 0 0\n\
ATOM 3 C C GLY A 1 2.8 0 0\nATOM 4 N N GLY A 2 4.1 0 0\n\
ATOM 5 C CA GLY A 2 5.5 0 0\nATOM 6 C C GLY A 2 6.9 0 0\n";

fn component() -> Component {
    let atom = |name: &str, element| ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: if name == "N" { -1 } else { 0 },
        aromatic: name == "N",
        leaving: false,
        stereo: (name == "CA").then_some(StereoConfiguration::S),
    };
    Component {
        id: "GLY".into(),
        name: "GLYCINE".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(b'G'),
        formula: Some("C2 H5 N O2".into()),
        atoms: vec![
            atom("N", Element::NITROGEN),
            atom("CA", Element::CARBON),
            atom("C", Element::CARBON),
        ]
        .into(),
        bonds: vec![
            ComponentBond {
                atom_a: "N".into(),
                atom_b: "CA".into(),
                order: BondOrder::Single,
                aromatic: false,
                stereo: None,
            },
            ComponentBond {
                atom_a: "CA".into(),
                atom_b: "C".into(),
                order: BondOrder::Single,
                aromatic: false,
                stereo: None,
            },
        ]
        .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn carbonyl_component() -> Component {
    let mut component = component();
    let mut atoms = component.atoms.to_vec();
    atoms[0].charge = 0;
    atoms[0].aromatic = false;
    for (name, element) in [("O", Element::OXYGEN), ("H", Element::HYDROGEN)] {
        atoms.push(ComponentAtom {
            name: name.into(),
            alternate_name: None,
            element,
            charge: 0,
            aromatic: false,
            leaving: false,
            stereo: None,
        });
    }
    component.atoms = atoms.into();
    let mut bonds = component.bonds.to_vec();
    for (first, second, order) in [("C", "O", BondOrder::Double), ("N", "H", BondOrder::Single)] {
        bonds.push(ComponentBond {
            atom_a: first.into(),
            atom_b: second.into(),
            order,
            aromatic: false,
            stereo: None,
        });
    }
    component.bonds = bonds.into();
    component
}

#[test]
fn linked_amide_nitrogens_stop_accepting_without_changing_terminal_roles_or_coordinates() {
    let (structure, _) = molframe_cif::read(
        &InputBuffer::from_bytes(STRUCTURE.as_bytes().to_vec()),
        &ReadOptions::new(),
    )
    .expect("fixture");
    let provider = MemoryProvider::new(DictionaryVersion::new("test"), [carbonyl_component()])
        .expect("dictionary");
    let policy = PolymerLinkPolicy::explicit(
        2.1,
        vec![PolymerLinkRule::new(
            ComponentKind::AminoAcid,
            "C",
            ComponentKind::AminoAcid,
            "N",
        )],
    );
    for (policy, second_accepts) in [(PolymerLinkPolicy::Disabled, true), (policy, false)] {
        let report = apply_component_chemistry(&structure, &provider, policy).expect("annotation");
        let annotation = |name| match report.structure.annotations().get(name) {
            Some(molframe_core::AtomAnnotation::Boolean(column)) => column,
            _ => panic!("boolean role column"),
        };
        assert_eq!(
            annotation(molframe_core::HBOND_ACCEPTOR_ANNOTATION).get(0),
            Some((true, molframe_core::Presence::Present))
        );
        assert_eq!(
            annotation(molframe_core::HBOND_ACCEPTOR_ANNOTATION).get(3),
            Some((second_accepts, molframe_core::Presence::Present))
        );
        assert_eq!(
            annotation(molframe_core::HBOND_DONOR_ANNOTATION).get(3),
            Some((true, molframe_core::Presence::Present))
        );
        assert_eq!(
            report.structure.positions().as_ptr(),
            structure.positions().as_ptr(),
            "chemistry never copies atom coordinates"
        );
    }
}

#[test]
fn deposited_covalent_links_refine_roles_even_when_polymer_inference_is_disabled() {
    let (structure, _) = molframe_cif::read(
        &InputBuffer::from_bytes(STRUCTURE.as_bytes().to_vec()),
        &ReadOptions::new(),
    )
    .expect("fixture");
    let mut data = structure.data().clone();
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: molframe_core::AtomIndex::new(2),
        atom_b: molframe_core::AtomIndex::new(3),
        order: BondOrder::Single,
        provenance: BondProvenance::File,
        kind: molframe_core::optional::OptionalSymbol::NONE,
    });
    data.bonds = bonds.finish();
    let structure = Structure::new(data);
    let provider = MemoryProvider::new(DictionaryVersion::new("test"), [carbonyl_component()])
        .expect("dictionary");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("annotation");
    let Some(molframe_core::AtomAnnotation::Boolean(column)) = report
        .structure
        .annotations()
        .get(molframe_core::HBOND_ACCEPTOR_ANNOTATION)
    else {
        panic!("acceptors");
    };
    assert_eq!(
        column.get(3),
        Some((false, molframe_core::Presence::Present))
    );
}

#[test]
fn unmatched_hydrogen_names_do_not_discard_matched_heavy_atom_chemistry() {
    let text = STRUCTURE.replace(
        "ATOM 4 N N GLY A 2",
        "ATOM 7 H H99 GLY A 1 0 1 0\nATOM 4 N N GLY A 2",
    );
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture");
    let provider =
        MemoryProvider::new(DictionaryVersion::new("test"), [component()]).expect("dictionary");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("annotation");
    let Some(molframe_core::AtomAnnotation::Integer(charges)) = report
        .structure
        .annotations()
        .get(molframe_core::FORMAL_CHARGE_ANNOTATION)
    else {
        panic!("charge column")
    };
    assert_eq!(
        charges.get(0),
        Some((-1, molframe_core::column::Presence::Present))
    );
    assert_eq!(
        charges.get(3),
        Some((0, molframe_core::column::Presence::Unknown))
    );
}

#[test]
fn ccd_edges_polymer_links_and_chain_classification_are_applied_together() {
    let input = InputBuffer::from_bytes(STRUCTURE.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider = MemoryProvider::new(DictionaryVersion::new("test-1"), [component()])
        .expect("component fixture is unique");
    let report = apply_component_chemistry(
        &structure,
        &provider,
        PolymerLinkPolicy::explicit(
            2.1,
            vec![PolymerLinkRule::new(
                ComponentKind::AminoAcid,
                "C",
                ComponentKind::AminoAcid,
                "N",
            )],
        ),
    )
    .expect("provider works");
    assert_eq!(report.structure.data().bonds.len(), 5);
    assert_eq!(
        report
            .structure
            .chain(molframe_core::ChainIndex::new(0))
            .map(molframe_core::structure::ChainRef::polymer_kind),
        Some(PolymerKind::Protein)
    );
    assert_eq!(report.dictionary_version.as_str(), "test-1");
    assert!(report.findings.is_empty());
    let Some(molframe_core::AtomAnnotation::Boolean(aromatic)) = report
        .structure
        .annotations()
        .get(molframe_core::AROMATIC_ATOM_ANNOTATION)
    else {
        panic!("aromatic annotation absent")
    };
    assert_eq!(
        aromatic.get(0),
        Some((true, molframe_core::Presence::Present))
    );
    let Some(molframe_core::AtomAnnotation::Integer(charges)) = report
        .structure
        .annotations()
        .get(molframe_core::FORMAL_CHARGE_ANNOTATION)
    else {
        panic!("formal charge annotation absent")
    };
    assert_eq!(charges.get(0), Some((-1, molframe_core::Presence::Present)));
}

#[test]
fn modelled_stereo_is_anchored_to_the_component_reference_geometry() {
    let source = "data_stereo\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
HETATM 1 C CTR LIG A 1 0 0 0\n\
HETATM 2 N A LIG A 1 1 0 0\n\
HETATM 3 O B LIG A 1 0 1 0\n\
HETATM 4 S C LIG A 1 0 0 -1\n";
    let atom = |name: &str, element, stereo| ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo,
    };
    let bond = |other: &str| ComponentBond {
        atom_a: "CTR".into(),
        atom_b: other.into(),
        order: BondOrder::Single,
        aromatic: false,
        stereo: None,
    };
    let component = Component {
        id: "LIG".into(),
        name: "CHIRAL LIGAND".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: vec![
            atom("CTR", Element::CARBON, Some(StereoConfiguration::R)),
            atom("A", Element::NITROGEN, None),
            atom("B", Element::OXYGEN, None),
            atom("C", Element::SULFUR, None),
        ]
        .into(),
        bonds: vec![bond("A"), bond("B"), bond("C")].into(),
        ideal_coordinates: Some(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
            .into(),
        ),
        model_coordinates: None,
    };
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider = MemoryProvider::new(DictionaryVersion::new("stereo"), [component])
        .expect("component fixture is unique");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("chemistry applies");
    let Some(molframe_core::AtomAnnotation::Symbol(stereo)) = report
        .structure
        .annotations()
        .get(molframe_core::STEREO_CONFIGURATION_ANNOTATION)
    else {
        panic!("stereo annotation absent")
    };
    let Some((configuration, molframe_core::Presence::Present)) = stereo.get(0) else {
        panic!("centre configuration absent")
    };
    assert_eq!(report.structure.resolve(configuration), Some("S"));
}

#[test]
fn an_unknown_component_is_reported_once_and_yields_a_known_empty_graph() {
    let input = InputBuffer::from_bytes(STRUCTURE.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider =
        MemoryProvider::new(DictionaryVersion::new("empty"), []).expect("empty fixture is unique");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("provider works");
    assert!(report.structure.data().bonds.is_available());
    assert_eq!(
        report
            .findings
            .iter()
            .filter(|finding| finding.code() == Code::W3201)
            .count(),
        1
    );
}

const HEADER: &str = "data_big\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

/// `chains` chains of `residues` glycines each; every residue lacks its carbon.
fn many_residues(chains: u32, residues: u32, component: &str) -> Structure {
    let mut text = String::from(HEADER);
    let mut id = 0;
    for chain in 0..chains {
        let label = char::from(b'A' + u8::try_from(chain).expect("few chains"));
        for seq in 1..=residues {
            for (name, element, x) in [("N", "N", 0.0_f32), ("CA", "C", 1.4)] {
                id += 1;
                writeln!(
                    text,
                    "ATOM {id} {element} {name} {component} {label} {seq} {x} {seq}.0 {chain}.0"
                )
                .expect("string write");
            }
        }
    }
    molframe_cif::read(
        &InputBuffer::from_bytes(text.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("fixture reads")
    .0
}

fn context<'a>(finding: &'a Diagnostic, label: &str) -> &'a str {
    finding
        .context()
        .iter()
        .find(|item| item.label() == label)
        .map_or("", |item| item.value())
}

fn fingerprint(report: &ChemistryReport) -> (Vec<String>, Vec<String>) {
    let bonds = report
        .structure
        .data()
        .bonds
        .iter()
        .map(|bond| format!("{bond:?}"))
        .collect();
    let findings = report.findings.iter().map(ToString::to_string).collect();
    (bonds, findings)
}

#[test]
fn repeated_missing_atom_findings_collapse_into_one_summary_with_exact_counts() {
    let structure = many_residues(2, 700, "GLY");
    let provider =
        MemoryProvider::new(DictionaryVersion::new("t"), [component()]).expect("dictionary");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("annotation");
    let missing: Vec<_> = report
        .findings
        .iter()
        .filter(|finding| finding.code() == Code::W3301)
        .collect();
    assert_eq!(missing.len(), 1, "one summary, not one line per residue");
    let summary = missing[0];
    assert_eq!(context(summary, "residues"), "1400");
    assert_eq!(context(summary, "components"), "1");
    assert_eq!(context(summary, "missing atoms"), "1400");
    assert_eq!(context(summary, "by component"), "GLYx1400");
    assert_eq!(context(summary, "sample").split(", ").count(), 8);
    assert!(context(summary, "sample").starts_with("GLY@0 (1 missing)"));
}

#[test]
fn inconsistent_residues_are_summarised_once_and_keep_their_severity() {
    let structure = many_residues(1, 300, "GLY");
    let mut odd = component();
    odd.atoms = odd
        .atoms
        .iter()
        .filter(|atom| atom.name.as_ref() != "CA")
        .cloned()
        .collect();
    odd.bonds = Arc::from([]);
    let provider = MemoryProvider::new(DictionaryVersion::new("t"), [odd]).expect("dictionary");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("annotation");
    let summary: Vec<_> = report
        .findings
        .iter()
        .filter(|finding| finding.code() == Code::W3202)
        .collect();
    assert_eq!(summary.len(), 1);
    assert_eq!(context(summary[0], "residues"), "300");
    assert_eq!(summary[0].severity(), Code::W3202.severity());
}

#[test]
fn parallel_units_match_a_single_threaded_pass_bond_for_bond() {
    // 600 residues per chain spans several work units per chain.
    let structure = many_residues(3, 600, "GLY");
    let provider =
        MemoryProvider::new(DictionaryVersion::new("t"), [component()]).expect("dictionary");
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool")
            .install(|| {
                apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
                    .expect("annotation")
            })
    };
    let sequential = run(1);
    let parallel = run(4);
    assert_eq!(fingerprint(&sequential), fingerprint(&parallel));
    assert_eq!(sequential.structure.data().bonds.len(), 3 * 600);
}

#[test]
fn a_lazy_dictionary_annotates_exactly_like_the_eager_one() {
    let ccd = "data_GLY\n_chem_comp.id GLY\n_chem_comp.name GLYCINE\n\
_chem_comp.type 'L-peptide linking'\n_chem_comp.one_letter_code G\n_chem_comp.formula 'C2 H5 N O2'\n\
loop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
_chem_comp_atom.pdbx_leaving_atom_flag\n_chem_comp_atom.pdbx_stereo_config\n\
N N 0 N N N\nCA C 0 N N S\nC C 0 N N N\n\
loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n\
_chem_comp_bond.pdbx_stereo_config\nN CA SING N N\nCA C SING N N\n\
data_UNUSED\n_chem_comp.id UNUSED\nthis is not a valid block\n";
    let structure = many_residues(1, 300, "GLY");
    let lazy = LazyCifProvider::new(
        InputBuffer::from_bytes(ccd.as_bytes().to_vec()),
        DictionaryVersion::new("t"),
    );
    assert_eq!(lazy.block_count(), 2);
    let eager =
        MemoryProvider::new(DictionaryVersion::new("t"), [component()]).expect("dictionary");
    let from_lazy = apply_component_chemistry(&structure, &lazy, PolymerLinkPolicy::Disabled)
        .expect("lazy annotation");
    let from_eager = apply_component_chemistry(&structure, &eager, PolymerLinkPolicy::Disabled)
        .expect("eager annotation");
    assert_eq!(lazy.lowered_count(), 1, "the unused block is never lowered");
    assert_eq!(fingerprint(&from_lazy).0, fingerprint(&from_eager).0);
}

#[path = "annotate/linkage_tests.rs"]
mod linkage_tests;
