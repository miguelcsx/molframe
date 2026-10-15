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
    let contact = contacts.row(0).expect("one contact");
    assert_eq!(contact.first.get(), 0);
    assert_eq!(contact.second.get(), 1);
    // C2-C3 are the closest atoms of the two residues, 0.5 Å apart.
    assert!((contact.min_distance - 0.5).abs() < 1e-5);
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

fn read(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn contact_pairs(structure: &Structure, min_separation: u32) -> Vec<(u32, u32)> {
    let Ok(map) = residue_contact_map(
        structure,
        2.0,
        min_separation,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    map.contacts()
        .iter()
        .map(|contact| (contact.first.get(), contact.second.get()))
        .collect()
}

// Residues 0 and 1 are adjacent in chain A; residue 2 is the first of chain B and
// sits 1 A from residue 1, so its global index is adjacent too. Residue 3 is the
// second of chain B, 1 A from residue 2.
const TWO_CHAINS: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 1 0 0\n\
ATOM 3 C CA GLY B 1 2 0 0\n\
ATOM 4 C CA GLY B 2 3 0 0\n";

#[test]
fn separation_one_excludes_adjacent_residues_of_the_same_chain() {
    let structure = read(TWO_CHAINS);
    // Residue indices differ by 1 within a chain: 0-1 and 2-3. Cross-chain 1-2
    // also differs by 1 globally but must survive. Pairs two apart (0-2, 1-3)
    // are 2 A away, within the cutoff.
    assert_eq!(
        contact_pairs(&structure, 0),
        vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    );
    assert_eq!(
        contact_pairs(&structure, 1),
        vec![(0, 2), (1, 2), (1, 3)],
        "same-chain neighbours are dropped and cross-chain pairs are kept"
    );
}

#[test]
fn cross_chain_pairs_survive_any_separation_filter() {
    let structure = read(TWO_CHAINS);
    assert_eq!(contact_pairs(&structure, 50), vec![(0, 2), (1, 2), (1, 3)]);
}
