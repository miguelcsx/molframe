use super::*;
use crate::{Component, ComponentAtom, ComponentBond, MemoryProvider, SnfgShape};
use molframe_core::contract::DictionaryVersion;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{AtomIndex, BondProvenance, BondRecord, Element};
use std::fmt::Write;
use std::sync::Arc;

fn component(size: usize, id: &str) -> Component {
    let names: Vec<_> = (0..size)
        .map(|i| {
            if i == 0 {
                "RO".to_owned()
            } else {
                format!("R{i}")
            }
        })
        .chain(["LEAVE".to_owned()])
        .collect();
    let atoms = names
        .iter()
        .enumerate()
        .map(|(i, name)| ComponentAtom {
            name: name.clone().into(),
            alternate_name: None,
            element: if i == 0 || i == size {
                Element::OXYGEN
            } else {
                Element::CARBON
            },
            charge: 0,
            aromatic: false,
            leaving: i == size,
            stereo: None,
        })
        .collect::<Vec<_>>();
    let bonds = (0..size)
        .map(|i| (i, (i + 1) % size))
        .chain([(1, size)])
        .map(|(a, b)| ComponentBond {
            atom_a: names[a].clone().into(),
            atom_b: names[b].clone().into(),
            order: BondOrder::Single,
            aromatic: false,
            stereo: None,
        })
        .collect::<Vec<_>>();
    Component {
        id: id.into(),
        name: "Test sugar".into(),
        kind: ComponentKind::Saccharide,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: Arc::from(atoms),
        bonds: Arc::from(bonds),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn fixture(
    size: usize,
    missing: Option<usize>,
    occupied: bool,
    terminal: Option<f32>,
    alternatives: bool,
) -> Structure {
    let mut text = String::from(
        "data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
",
    );
    let mut count = 0;
    for alt in if alternatives {
        vec!["A", "B"]
    } else {
        vec!["."]
    } {
        for i in 0..size {
            if missing == Some(i) || (alternatives && i == 0 && alt == "B") {
                continue;
            }
            let name = if i == 0 {
                "RO".to_owned()
            } else {
                format!("R{i}")
            };
            let angle = (f32::from(u8::try_from(i).unwrap()) - 1.0) * std::f32::consts::TAU
                / f32::from(u8::try_from(size).unwrap());
            count += 1;
            let alt = if i == 0 { "." } else { alt };
            writeln!(
                text,
                "HETATM {count} {} {name} {alt} ZZZ A 1 {} {} 0",
                if i == 0 { "O" } else { "C" },
                1.5 * angle.cos(),
                1.5 * angle.sin()
            )
            .unwrap();
        }
    }
    if occupied {
        count += 1;
        writeln!(text, "HETATM {count} O LEAVE . ZZZ A 1 2.8 0 0").unwrap();
    }
    if let Some(distance) = terminal {
        count += 1;
        writeln!(text, "ATOM {count} N ND2 . ASN B 1 {} 0 0", 1.5 + distance).unwrap();
        count += 1;
        writeln!(text, "ATOM {count} C CG . ASN B 1 {} 0 0", 2.8 + distance).unwrap();
    }
    let (structure, _) = molframe_cif::read(
        &InputBuffer::from_bytes(text.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("synthetic CIF");
    if terminal.is_none() {
        return structure;
    }
    let mut data = structure.data().clone();
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: AtomIndex::new(count - 2),
        atom_b: AtomIndex::new(count - 1),
        order: BondOrder::Single,
        provenance: BondProvenance::File,
    });
    data.bonds = bonds.finish();
    Structure::new(data)
}

fn perceive(structure: &Structure, size: usize) -> CarbohydrateReport {
    let provider = MemoryProvider::new(DictionaryVersion::new("test"), [component(size, "ZZZ")])
        .expect("provider");
    carbohydrates(structure, Some(&provider), CarbohydrateOptions::default()).expect("inventory")
}

#[test]
fn curated_symbols_preserve_molstar_shapes_colors_and_aliases() {
    assert_eq!(snfg_symbol("NAG").unwrap().abbreviation, "GlcNAc");
    assert_eq!(snfg_symbol("MAN"), snfg_symbol("BMA"));
    assert_eq!(snfg_symbol("FUC").unwrap().color, 0xed_1c_24);
    assert_eq!(
        snfg_symbol("GalN").unwrap().secondary_color,
        Some(0xf1_ec_e1)
    );
    assert_eq!(snfg_symbol("GCU").unwrap().shape, SnfgShape::DividedDiamond);
    assert!(snfg_symbol("0GA").is_none());
    assert!(snfg_symbol("ZZZ").is_none());
}

#[test]
fn ring_topology_identifies_both_furanose_and_pyranose_without_name_templates() {
    for size in [5, 6] {
        let report = perceive(&fixture(size, None, false, None, false), size);
        assert_eq!(report.monosaccharides.len(), 1);
        let sugar = &report.monosaccharides[0];
        assert_eq!(sugar.ring_atoms.len(), size);
        assert_eq!(sugar.symbol.abbreviation, "Unk");
        assert_eq!(sugar.anomeric_atom, Some(AtomIndex::new(1)));
        let geometry = sugar.geometry.expect("nondegenerate ring");
        assert!((geometry.normal[2].abs() - 1.0).abs() < 1e-6);
        assert_eq!(report.dictionary_version.as_ref().unwrap().as_str(), "test");
    }
}

#[test]
fn missing_ring_atoms_are_reported_and_proximity_never_fabricates_a_ring() {
    let source = fixture(6, Some(3), false, None, false);
    let report = perceive(&source, 6);
    assert!(report.monosaccharides.is_empty());
    assert_eq!(report.incomplete_residues.len(), 1);
    let source = fixture(6, None, false, None, false);
    assert!(
        carbohydrates(&source, None, CarbohydrateOptions::default())
            .unwrap()
            .monosaccharides
            .is_empty()
    );
}

#[test]
fn shared_atoms_do_not_mix_alternate_conformers() {
    let source = fixture(6, None, false, None, true);
    let report = perceive(&source, 6);
    assert_eq!(report.monosaccharides.len(), 2);
    for sugar in report.monosaccharides {
        let labels: std::collections::BTreeSet<_> = sugar
            .ring_atoms
            .iter()
            .filter_map(|atom| source.atom(*atom).unwrap().alt_label())
            .collect();
        assert_eq!(labels.len(), 1);
    }
}

#[test]
fn fallback_requires_vacant_anomeric_chemistry_and_obeys_the_two_angstrom_boundary() {
    for (distance, expected) in [(1.4, 1), (2.0, 1), (2.001, 0), (0.8, 0)] {
        let source = fixture(6, None, false, Some(distance), false);
        let report = perceive(&source, 6);
        assert_eq!(report.terminal_links.len(), expected, "distance {distance}");
        if expected == 1 {
            assert_eq!(
                report.terminal_links[0].provenance,
                BondProvenance::InferredDistance
            );
        }
    }
    assert!(
        perceive(&fixture(6, None, true, Some(1.4), false), 6)
            .terminal_links
            .is_empty()
    );
}

#[test]
fn declared_attachment_keeps_file_provenance_and_source_snapshot_unchanged() {
    let source = fixture(6, None, false, Some(1.4), false);
    let mut data = source.data().clone();
    let mut bonds = BondTableBuilder::new();
    for bond in data.bonds.iter() {
        bonds.push(bond);
    }
    bonds.push(BondRecord {
        atom_a: AtomIndex::new(1),
        atom_b: AtomIndex::new(6),
        order: BondOrder::Single,
        provenance: BondProvenance::File,
    });
    data.bonds = bonds.finish();
    let source = Structure::new(data);
    let report = perceive(&source, 6);
    assert_eq!(report.terminal_links.len(), 1);
    assert_eq!(report.terminal_links[0].provenance, BondProvenance::File);
    assert_eq!(source.data().bonds.len(), 2);
}

#[test]
fn competing_donors_and_saturated_or_unconnected_acceptors_are_refused() {
    let source = fixture(6, None, false, Some(1.4), true);
    let report = perceive(&source, 6);
    assert_eq!(report.monosaccharides.len(), 2);
    assert!(report.terminal_links.is_empty());
    let source = fixture(6, None, false, Some(1.4), false);
    let mut data = source.data().clone();
    data.bonds = BondTableBuilder::new().finish();
    assert!(perceive(&Structure::new(data), 6).terminal_links.is_empty());
    let mut data = source.data().clone();
    let mut bonds = BondTableBuilder::new();
    for bond in data.bonds.iter() {
        bonds.push(bond);
    }
    bonds.push(BondRecord {
        atom_a: AtomIndex::new(6),
        atom_b: AtomIndex::new(3),
        order: BondOrder::Single,
        provenance: BondProvenance::File,
    });
    data.bonds = bonds.finish();
    assert!(perceive(&Structure::new(data), 6).terminal_links.is_empty());
}

#[test]
fn existing_inferred_edges_obey_vacancy_and_distance_instead_of_becoming_declarations() {
    for (occupied, distance, expected) in [(false, 1.4, 1), (true, 1.4, 0), (false, 2.01, 0)] {
        let source = fixture(6, None, occupied, Some(distance), false);
        let mut data = source.data().clone();
        let mut bonds = BondTableBuilder::new();
        for bond in data.bonds.iter() {
            bonds.push(bond);
        }
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(1),
            atom_b: AtomIndex::new(if occupied { 7 } else { 6 }),
            order: BondOrder::Single,
            provenance: BondProvenance::InferredDistance,
        });
        data.bonds = bonds.finish();
        let report = perceive(&Structure::new(data), 6);
        assert_eq!(report.terminal_links.len(), expected);
        assert!(
            report
                .terminal_links
                .iter()
                .all(|link| link.provenance == BondProvenance::InferredDistance)
        );
    }
}

#[path = "real_tests.rs"]
mod real;
