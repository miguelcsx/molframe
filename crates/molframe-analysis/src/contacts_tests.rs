use super::{atom_contacts, visit_atom_contacts};
use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;

// Two residues on the x axis: C1(0) C2(1) in residue 1, C3(1.5) C4(10) in
// residue 2. The only separations under 2 Å are C1-C2, C1-C3 and C2-C3.
const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 1 0 0\n\
ATOM 3 C C3 LIG A 2 1.5 0 0\n\
ATOM 4 C C4 LIG A 2 10 0 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn a_cutoff_reports_exactly_the_near_pairs_in_order() {
    let structure = structure();
    let Ok(contacts) = atom_contacts(
        &structure,
        2.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid cutoff");
    };
    let pairs: Vec<(u32, u32)> = contacts
        .iter()
        .map(|contact| (contact.first.get(), contact.second.get()))
        .collect();
    assert_eq!(pairs, vec![(0, 1), (0, 2), (1, 2)]);
}

#[test]
fn the_reported_distance_matches_the_geometry() {
    let structure = structure();
    let Ok(contacts) = atom_contacts(
        &structure,
        2.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid cutoff");
    };
    // C2-C3 are 0.5 Å apart.
    let touching = contacts
        .iter()
        .find(|contact| contact.first.get() == 1 && contact.second.get() == 2);
    let Some(touching) = touching else {
        panic!("expected the C2-C3 contact");
    };
    assert!((touching.distance - 0.5).abs() < 1e-5);
}

#[test]
fn every_backend_agrees_on_the_contacts() {
    let structure = structure();
    let baseline = atom_contacts(
        &structure,
        2.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    );
    let Ok(baseline) = baseline else {
        panic!("valid");
    };
    for backend in [SpatialBackend::CellList, SpatialBackend::KdTree] {
        let Ok(other) = atom_contacts(&structure, 2.0, backend, &ExecutionContext::default())
        else {
            panic!("valid");
        };
        assert_eq!(baseline, other, "backend {backend:?} disagreed");
    }
}

#[test]
fn a_negative_cutoff_is_rejected() {
    let structure = structure();
    assert!(
        atom_contacts(
            &structure,
            -1.0,
            SpatialBackend::BruteForce,
            &ExecutionContext::default(),
        )
        .is_err()
    );
}

#[test]
fn visitor_emits_the_materialized_contact_set_without_retaining_it() {
    let structure = structure();
    let context = ExecutionContext::default();
    let expected = match atom_contacts(&structure, 2.0, SpatialBackend::CellList, &context) {
        Ok(contacts) => contacts,
        Err(error) => panic!("materialized contacts failed: {error}"),
    };
    let mut observed = Vec::new();
    let visited = visit_atom_contacts(
        &structure,
        2.0,
        SpatialBackend::CellList,
        &context,
        |contact| observed.push(contact),
    );
    assert!(visited.is_ok());
    observed.sort_unstable_by_key(|contact| (contact.first, contact.second));
    assert_eq!(observed, expected);
}

#[test]
fn worker_count_does_not_change_the_contacts() {
    let structure = structure();
    let serial = match atom_contacts(
        &structure,
        2.0,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    ) {
        Ok(contacts) => contacts,
        Err(error) => panic!("serial contacts failed: {error}"),
    };
    for workers in [1, 2, 4, 8] {
        let context = match ExecutionContext::builder().worker_budget(workers).build() {
            Ok(context) => context,
            Err(error) => panic!("valid execution context: {error}"),
        };
        let parallel = match atom_contacts(&structure, 2.0, SpatialBackend::CellList, &context) {
            Ok(contacts) => contacts,
            Err(error) => panic!("parallel contacts failed: {error}"),
        };
        assert_eq!(
            serial, parallel,
            "worker count {workers} changed the contacts"
        );
    }
}
