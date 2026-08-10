use super::{PolymerRoleProfile, PolymerRoleRule, apply_polymer_role_profile};
use crate::{Component, ComponentAtom, ComponentKind, MemoryProvider, PolymerAtomRole};
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::{AtomAnnotation, Element, Presence};
use std::sync::Arc;

#[test]
fn caller_profile_assigns_roles_without_a_library_name_table() {
    let input = InputBuffer::from_bytes(
        b"data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 N Q GLY A 1 0 0 0\n".to_vec(),
    );
    let Ok((structure, _)) = pdbiox_cif::read(&input, &ReadOptions::new()) else {
        panic!("valid structure");
    };
    let provider = MemoryProvider::new(DictionaryVersion::new("ccd-test"), [component()])
        .expect("component fixture is unique");
    let profile = PolymerRoleProfile {
        id: "roles-test-v1".into(),
        rules: Arc::from([PolymerRoleRule {
            component_id: None,
            component_kind: Some(ComponentKind::AminoAcid),
            atom_name: "Q".into(),
            role: PolymerAtomRole::PROTEIN_NITROGEN,
        }]),
    };
    let Ok(report) = apply_polymer_role_profile(&structure, &provider, &profile) else {
        panic!("valid explicit role profile");
    };
    let Some(AtomAnnotation::Integer(roles)) = report
        .structure
        .annotations()
        .get(pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION)
    else {
        panic!("role annotation should be present");
    };
    assert_eq!(
        roles.get(0),
        Some((PolymerAtomRole::PROTEIN_NITROGEN.code(), Presence::Present))
    );
}

fn component() -> Component {
    Component {
        id: "GLY".into(),
        name: "test".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(b'G'),
        formula: None,
        atoms: Arc::from([ComponentAtom {
            name: "Q".into(),
            alternate_name: None,
            element: Element::NITROGEN,
            charge: 0,
            aromatic: false,
            leaving: false,
            stereo: None,
        }]),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}
