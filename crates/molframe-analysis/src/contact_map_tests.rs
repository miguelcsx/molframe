use super::residue_contact_map;
use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;

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
fn adjacent_residues_touch_at_their_closest_atoms() {
    let structure = structure();
    let Ok(map) = residue_contact_map(
        &structure,
        2.0,
        0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert_eq!(map.residue_count(), 2);
    let contacts = map.contacts();
    assert_eq!(contacts.len(), 1);
    assert_eq!(contacts[0].first.get(), 0);
    assert_eq!(contacts[0].second.get(), 1);
    // C2-C3 are the closest atoms of the two residues, 0.5 Å apart.
    assert!((contacts[0].min_distance - 0.5).abs() < 1e-5);
}

#[test]
fn a_separation_filter_drops_neighbouring_residues() {
    let structure = structure();
    let Ok(map) = residue_contact_map(
        &structure,
        2.0,
        2,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(
        map.contacts().is_empty(),
        "the only pair is one residue apart"
    );
}

#[test]
fn intra_residue_contacts_are_never_reported() {
    // A cutoff that also captures the C1-C2 intra-residue pair must still yield
    // just the one inter-residue contact.
    let structure = structure();
    let Ok(map) = residue_contact_map(
        &structure,
        2.0,
        0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(
        map.contacts().iter().all(|c| c.first != c.second),
        "a residue cannot contact itself"
    );
}

#[test]
fn worker_count_does_not_change_the_map() {
    let structure = structure();
    let serial = match residue_contact_map(
        &structure,
        2.0,
        0,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    ) {
        Ok(map) => map,
        Err(error) => panic!("serial map failed: {error}"),
    };
    for workers in [1, 2, 4, 8] {
        let context = match ExecutionContext::builder().worker_budget(workers).build() {
            Ok(context) => context,
            Err(error) => panic!("valid execution context: {error}"),
        };
        let parallel =
            match residue_contact_map(&structure, 2.0, 0, SpatialBackend::CellList, &context) {
                Ok(map) => map,
                Err(error) => panic!("parallel map failed: {error}"),
            };
        assert_eq!(serial, parallel, "worker count {workers} changed the map");
    }
}
