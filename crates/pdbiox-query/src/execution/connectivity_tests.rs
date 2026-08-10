use super::*;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};

const SOURCE: &str = "data_bonds\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C A GLY A 1 0 0 0\nATOM 2 C B GLY A 1 1 0 0\n\
ATOM 3 C C GLY A 1 2 0 0\nATOM 4 C D GLY A 1 3 0 0\n\
ATOM 5 C E GLY A 1 4 0 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = pdbiox_cif::read(&input, &ReadOptions::new()).expect("fixture parses");
    let mut bonds = BondTableBuilder::new();
    for (a, b) in [(0, 1), (1, 2), (3, 4)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(a),
            atom_b: AtomIndex::new(b),
            order: BondOrder::Single,
            provenance: BondProvenance::File,
        });
    }
    let mut data = structure.data().clone();
    data.bonds = bonds.finish();
    Structure::new(data)
}

#[test]
fn bounded_bond_traversal_includes_each_depth_and_respects_the_universe() {
    let structure = structure();
    let target = AtomSelection::from_sorted(vec![0]);
    let universe = AtomSelection::range(1..5);
    assert_eq!(
        bonded(&structure, &universe, &target, 2).expect("graph exists"),
        AtomSelection::range(1..3)
    );
}

#[test]
fn fragment_expansion_selects_only_components_touched_by_the_target() {
    let structure = structure();
    let target = AtomSelection::from_sorted(vec![1]);
    assert_eq!(
        same_fragment(
            &structure,
            &AtomSelection::All(structure.atom_count()),
            &target
        )
        .expect("graph exists"),
        AtomSelection::range(0..3)
    );
}
