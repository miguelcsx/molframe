use super::{WaterBridgeOptions, water_bridges};
use crate::HydrogenBondOptions;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::{
    AnnotationColumn, AtomAnnotation, BondOrder, BondProvenance, BondRecord, BondTableBuilder,
    ExecutionContext, Presence,
};
use pdbiox_spatial::SpatialBackend;

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N DON A 1 0 0 0\nATOM 2 H H DON A 1 1 0 0\n\
HETATM 3 O O HOH W 2 2.8 0 0\nHETATM 4 H H1 HOH W 2 3.8 0 0\n\
ATOM 5 O O ACC B 3 5.6 0 0\n";

#[test]
fn two_oriented_hydrogen_bonds_form_one_water_bridge() {
    let bridges = water_bridges(
        &annotated_structure(),
        WaterBridgeOptions {
            hydrogen_bonds: HydrogenBondOptions {
                maximum_donor_acceptor_distance: 3.5,
                minimum_angle_degrees: 150.0,
                backend: SpatialBackend::BruteForce,
                periodic: false,
            },
        },
        &ExecutionContext::default(),
    )
    .expect("valid bridge network");
    assert_eq!(bridges.len(), 1);
    assert_eq!(bridges[0].water.get(), 2);
    assert_eq!(bridges[0].first.get(), 0);
    assert_eq!(bridges[0].second.get(), 4);
}

fn annotated_structure() -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let structure = match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let mut data = structure.data().clone();
    data.annotations.insert(
        pdbiox_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(boolean_roles(5, &[0, 2])),
    );
    data.annotations.insert(
        pdbiox_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(boolean_roles(5, &[2, 4])),
    );
    data.annotations.insert(
        pdbiox_core::COMPONENT_KIND_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_entries((0..5).map(|atom| {
                let kind = if [2, 3].contains(&atom) {
                    pdbiox_chem::ComponentKind::Solvent
                } else {
                    pdbiox_chem::ComponentKind::AminoAcid
                };
                (kind.code(), Presence::Present)
            }))
            .expect("small annotation column"),
        ),
    );
    let mut bonds = BondTableBuilder::new();
    for (first, second) in [(0, 1), (2, 3)] {
        bonds.push(BondRecord {
            atom_a: pdbiox_core::AtomIndex::new(first),
            atom_b: pdbiox_core::AtomIndex::new(second),
            order: BondOrder::Single,
            provenance: BondProvenance::ChemicalComponentDictionary,
        });
    }
    data.bonds = bonds.finish();
    pdbiox_core::Structure::new(data)
}

fn boolean_roles(atoms: u32, selected: &[u32]) -> AnnotationColumn<bool> {
    AnnotationColumn::from_entries((0..atoms).map(|atom| {
        if selected.contains(&atom) {
            (true, Presence::Present)
        } else {
            (false, Presence::Inapplicable)
        }
    }))
    .expect("small annotation column")
}
