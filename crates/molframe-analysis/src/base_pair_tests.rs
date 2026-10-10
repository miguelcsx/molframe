use super::{BasePairOptions, WatsonCrickGeometry, base_pairs};
use crate::HydrogenBondOptions;
use molframe_chem::{
    Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider, PolymerLinkPolicy,
    apply_component_chemistry,
};
use molframe_core::contract::DictionaryVersion;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{BondOrder, Element, ExecutionContext};
use molframe_spatial::SpatialBackend;
use std::fmt::Write as _;
use std::sync::Arc;

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N6 ADE A 1 0 0 0\n\
ATOM 2 H H61 ADE A 1 1 0 0\n\
ATOM 3 O O4 URA B 1 2.8 0 0\n";

#[test]
fn ccd_identity_and_oriented_bond_define_a_pair() {
    let structure = annotated_structure();
    let Ok(pairs) = base_pairs(
        &structure,
        &provider(b'U'),
        options(),
        &ExecutionContext::default(),
    ) else {
        panic!("valid base-pair analysis");
    };
    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0].hydrogen_bond_count, 1);
}

#[test]
fn non_complementary_ccd_codes_do_not_pair() {
    let structure = annotated_structure();
    let Ok(pairs) = base_pairs(
        &structure,
        &provider(b'G'),
        options(),
        &ExecutionContext::default(),
    ) else {
        panic!("valid non-complementary analysis");
    };
    assert!(pairs.is_empty());
}

#[test]
fn supporting_bond_count_is_an_explicit_policy() {
    let structure = annotated_structure();
    let mut options = options();
    options.minimum_hydrogen_bonds = 2;
    let Ok(pairs) = base_pairs(
        &structure,
        &provider(b'U'),
        options,
        &ExecutionContext::default(),
    ) else {
        panic!("valid stringent analysis");
    };
    assert!(pairs.is_empty());
}

fn annotated_structure() -> molframe_core::Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    match apply_component_chemistry(&structure, &provider(b'U'), PolymerLinkPolicy::Disabled) {
        Ok(report) => report.structure,
        Err(error) => panic!("chemistry failed: {error}"),
    }
}

fn options() -> BasePairOptions {
    BasePairOptions {
        hydrogen_bonds: HydrogenBondOptions {
            maximum_donor_acceptor_distance: 3.5,
            minimum_angle_degrees: 150.0,
            backend: SpatialBackend::BruteForce,
            periodic: false,
        },
        minimum_hydrogen_bonds: 1,
        geometry: None,
    }
}

fn provider(second_code: u8) -> MemoryProvider {
    MemoryProvider::new(
        DictionaryVersion::new("test"),
        [donor_component(), acceptor_component(second_code)],
    )
    .expect("component fixtures are unique")
}

fn donor_component() -> Component {
    Component {
        id: "ADE".into(),
        name: "adenine test".into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(b'A'),
        formula: None,
        atoms: Arc::from([
            atom("N6", Element::NITROGEN),
            atom("H61", Element::HYDROGEN),
        ]),
        bonds: Arc::from([ComponentBond {
            atom_a: "N6".into(),
            atom_b: "H61".into(),
            order: BondOrder::Single,
            aromatic: false,
            stereo: None,
        }]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn acceptor_component(code: u8) -> Component {
    Component {
        id: "URA".into(),
        name: "acceptor test".into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(code),
        formula: None,
        atoms: Arc::from([atom("O4", Element::OXYGEN)]),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn atom(name: &str, element: Element) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    }
}

#[test]
fn a_contact_off_the_watson_crick_edge_is_not_a_pair() {
    // Adenine N7 (Hoogsteen edge) accepting from uracil N3-H3.
    let hoogsteen = geometry_fixture(&[
        ("ADE", "N7", Element::NITROGEN, [5.0, 0.0, 0.0]),
        ("URA", "N3", Element::NITROGEN, [8.0, 0.0, 0.0]),
        ("URA", "H3", Element::HYDROGEN, [7.0, 0.0, 0.0]),
    ]);
    assert_eq!(bond_count(&hoogsteen), 1, "the contact is a real H-bond");
    assert!(pairs_of(&hoogsteen, Some(WatsonCrickGeometry::default())).is_empty());
    assert!(pairs_of(&hoogsteen, None).is_empty());
}

#[test]
fn a_planar_edge_bond_at_pair_width_is_accepted() {
    let canonical = geometry_fixture(&[
        ("ADE", "N6", Element::NITROGEN, [5.0, 0.0, 0.0]),
        ("ADE", "H61", Element::HYDROGEN, [6.0, 0.0, 0.0]),
        ("URA", "O4", Element::OXYGEN, [8.0, 0.0, 0.0]),
    ]);
    assert_eq!(bond_count(&canonical), 1);
    assert_eq!(pairs_of(&canonical, None).len(), 1, "edge test alone");
    let pairs = pairs_of(&canonical, Some(WatsonCrickGeometry::default()));
    assert_eq!(pairs.len(), 1);
}

#[test]
fn a_pair_with_wrong_c1_distance_or_tilted_planes_is_rejected() {
    let edge = [
        ("ADE", "N6", Element::NITROGEN, [5.0, 0.0, 0.0]),
        ("ADE", "H61", Element::HYDROGEN, [6.0, 0.0, 0.0]),
        ("URA", "O4", Element::OXYGEN, [8.0, 0.0, 0.0]),
    ];
    let limits = WatsonCrickGeometry::default();
    let narrow = WatsonCrickGeometry {
        minimum_c1_distance: 11.0,
        ..limits
    };
    assert!(pairs_of(&geometry_fixture(&edge), Some(narrow)).is_empty());
    let strict_planes = WatsonCrickGeometry {
        maximum_plane_angle_degrees: 0.0,
        ..limits
    };
    // Exactly coplanar rings still pass a zero-degree limit.
    assert_eq!(
        pairs_of(&geometry_fixture(&edge), Some(strict_planes)).len(),
        1
    );
    let tilted = geometry_fixture_with(&edge, 60.0);
    assert!(pairs_of(&tilted, Some(limits)).is_empty());
}

type Fixture = (molframe_core::Structure, MemoryProvider);

fn pairs_of(fixture: &Fixture, geometry: Option<WatsonCrickGeometry>) -> Vec<super::BasePair> {
    let mut options = options();
    options.geometry = geometry;
    match base_pairs(
        &fixture.0,
        &fixture.1,
        options,
        &ExecutionContext::default(),
    ) {
        Ok(pairs) => pairs,
        Err(error) => panic!("analysis failed: {error}"),
    }
}

fn bond_count(fixture: &Fixture) -> usize {
    match crate::hydrogen_bonds(
        &fixture.0,
        options().hydrogen_bonds,
        &ExecutionContext::default(),
    ) {
        Ok(bonds) => bonds.iter().count(),
        Err(error) => panic!("hydrogen bonds failed: {error}"),
    }
}

fn geometry_fixture(extra: &[(&str, &str, Element, [f32; 3])]) -> Fixture {
    geometry_fixture_with(extra, 0.0)
}

/// Two coplanar ring triangles (C2, C4, C6) 10.5 Å apart at the C1' atoms, with
/// the uracil rotated by `tilt_degrees` about the x axis, plus `extra` atoms.
fn geometry_fixture_with(extra: &[(&str, &str, Element, [f32; 3])], tilt_degrees: f32) -> Fixture {
    let (sin, cos) = tilt_degrees.to_radians().sin_cos();
    let tilt = |p: [f32; 3]| [p[0], p[1] * cos - p[2] * sin, p[1] * sin + p[2] * cos];
    let mut atoms: Vec<(&str, &str, Element, [f32; 3])> = vec![
        ("ADE", "C1'", Element::CARBON, [0.0, 0.0, 0.0]),
        ("ADE", "C2", Element::CARBON, [2.0, 1.0, 0.0]),
        ("ADE", "C4", Element::CARBON, [3.0, -1.0, 0.0]),
        ("ADE", "C6", Element::CARBON, [4.0, 1.0, 0.0]),
        ("URA", "C1'", Element::CARBON, [10.5, 0.0, 0.0]),
        ("URA", "C2", Element::CARBON, tilt([9.0, 1.0, 0.0])),
        ("URA", "C4", Element::CARBON, tilt([8.0, -1.0, 0.0])),
        ("URA", "C6", Element::CARBON, tilt([7.0, 1.0, 0.0])),
    ];
    atoms.extend_from_slice(extra);
    // Residues must be contiguous rows.
    atoms.sort_by_key(|(comp, ..)| *comp);
    let mut source = String::from(
        "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    for (id, (comp, name, element, p)) in atoms.iter().enumerate() {
        let chain = if *comp == "ADE" { "A" } else { "B" };
        let _ = writeln!(
            source,
            "ATOM {} {} {name} {comp} {chain} 1 {} {} {}",
            id + 1,
            element.symbol(),
            p[0],
            p[1],
            p[2]
        );
    }
    let component = |id: &str, code: u8| Component {
        id: id.into(),
        name: id.into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(code),
        formula: None,
        atoms: atoms
            .iter()
            .filter(|(comp, ..)| *comp == id)
            .map(|(_, name, element, _)| atom(name, *element))
            .collect::<Vec<_>>()
            .into(),
        bonds: donor_bonds(id, &atoms),
        ideal_coordinates: None,
        model_coordinates: None,
    };
    let provider = MemoryProvider::new(
        DictionaryVersion::new("test"),
        [component("ADE", b'A'), component("URA", b'U')],
    )
    .expect("component fixtures are unique");
    let input = InputBuffer::from_bytes(source.into_bytes());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let structure =
        match apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled) {
            Ok(report) => report.structure,
            Err(error) => panic!("chemistry failed: {error}"),
        };
    (structure, provider)
}

fn donor_bonds(id: &str, atoms: &[(&str, &str, Element, [f32; 3])]) -> Arc<[ComponentBond]> {
    let mut bonds = Vec::new();
    for (donor, hydrogen) in [("N6", "H61"), ("N3", "H3")] {
        let has = |wanted: &str| atoms.iter().any(|(c, n, ..)| *c == id && *n == wanted);
        if has(donor) && has(hydrogen) {
            bonds.push(ComponentBond {
                atom_a: donor.into(),
                atom_b: hydrogen.into(),
                order: BondOrder::Single,
                aromatic: false,
                stereo: None,
            });
        }
    }
    bonds.into()
}
