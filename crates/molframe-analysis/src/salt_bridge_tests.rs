use super::{SaltBridgeError, SaltBridgeOptions, salt_bridges, salt_bridges_with_options};
use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;
use std::fmt::Write;

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

#[test]
fn an_aspartate_oxygen_near_a_lysine_nitrogen_is_a_bridge() {
    let source = format!(
        "{HEADER}\
ATOM 1 O OD1 ASP A 1 0 0 0\n\
ATOM 2 N NZ LYS A 2 3.5 0 0\n"
    );
    let structure = crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, 1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert_eq!(bridges.len(), 1);
    let bridge = bridges.row(0).expect("one bridge");
    assert_eq!(bridge.anion.get(), 0);
    assert_eq!(bridge.cation.get(), 1);
    assert!((bridge.distance - 3.5).abs() < 1e-5);
}

#[test]
fn a_pair_beyond_the_cutoff_is_not_a_bridge() {
    let source = format!(
        "{HEADER}\
ATOM 1 O OD1 ASP A 1 0 0 0\n\
ATOM 2 N NZ LYS A 2 5 0 0\n"
    );
    let structure = crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, 1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(bridges.is_empty());
}

#[test]
fn two_anions_do_not_bridge_each_other() {
    let source = format!(
        "{HEADER}\
ATOM 1 O OD1 ASP A 1 0 0 0\n\
ATOM 2 O OE1 GLU A 2 3 0 0\n"
    );
    let structure =
        crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, -1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(bridges.is_empty(), "like charges do not bridge");
}

#[test]
fn a_non_charged_side_chain_atom_is_ignored() {
    // Aspartate's CB carbon is not part of the carboxylate.
    let source = format!(
        "{HEADER}\
ATOM 1 C CB ASP A 1 0 0 0\n\
ATOM 2 N NZ LYS A 2 3 0 0\n"
    );
    let structure = crate::chemistry_test_support::charges(&structure(&source), &[(1, 1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(bridges.is_empty());
}

#[test]
fn the_bridge_list_is_identical_at_every_worker_count() {
    // Alternating charged atoms across a lattice, dense enough that the search
    // spans several query blocks.
    let mut source = String::from(HEADER);
    let mut charges = Vec::new();
    let mut serial_id = 0u32;
    for x in 0..7i16 {
        for y in 0..7i16 {
            for z in 0..7i16 {
                let anion = serial_id.is_multiple_of(2);
                let (name, element) = if anion { ("OD1", "O") } else { ("NZ", "N") };
                serial_id += 1;
                writeln!(
                    source,
                    "ATOM {serial_id} {element} {name} LIG A {serial_id} {:.3} {:.3} {:.3}",
                    f32::from(x) * 3.0,
                    f32::from(y) * 3.0,
                    f32::from(z) * 3.0
                )
                .expect("fixture row");
                charges.push((serial_id - 1, if anion { -1 } else { 1 }));
            }
        }
    }

    let structure = crate::chemistry_test_support::charges(&structure(&source), &charges);
    let Ok(serial) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    ) else {
        panic!("valid lattice");
    };
    assert!(!serial.is_empty(), "the fixture must produce bridges");

    for workers in [2, 4, 8, 16] {
        let context = match ExecutionContext::builder().worker_budget(workers).build() {
            Ok(context) => context,
            Err(error) => panic!("valid execution context: {error}"),
        };
        let Ok(parallel) = salt_bridges(&structure, 4.0, SpatialBackend::CellList, &context) else {
            panic!("valid lattice");
        };
        assert_eq!(
            parallel, serial,
            "worker count {workers} changed the result"
        );
    }
}

fn pair(cell: &str, residues: (u32, u32), second_x: f32) -> Structure {
    let header = HEADER.replacen("data_s\n", &format!("data_s\n{cell}"), 1);
    let source = format!(
        "{header}\
ATOM 1 O OD1 ASP A {} 5 5 5\n\
ATOM 2 N NZ LYS A {} {second_x} 5 5\n",
        residues.0, residues.1
    );
    crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, 1)])
}

fn count(structure: &Structure, options: SaltBridgeOptions) -> Result<usize, SaltBridgeError> {
    salt_bridges_with_options(structure, options, &ExecutionContext::default())
        .map(|table| table.len())
}

#[test]
fn same_residue_pairs_are_excluded_unless_included() {
    let same = pair("", (1, 1), 8.5);
    let mut options = SaltBridgeOptions::new(4.0, SpatialBackend::BruteForce);
    assert_eq!(count(&same, options), Ok(0));
    options.include_same_residue_pairs = true;
    assert_eq!(count(&same, options), Ok(1));
    assert_eq!(
        salt_bridges(
            &same,
            4.0,
            SpatialBackend::BruteForce,
            &ExecutionContext::default()
        )
        .map(|t| t.len()),
        Ok(0)
    );
}

#[test]
fn bonded_pairs_are_excluded_unless_included() {
    use molframe_core::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
    let base = pair("", (1, 2), 6.5);
    let mut data = base.data().clone();
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: molframe_core::AtomIndex::new(0),
        atom_b: molframe_core::AtomIndex::new(1),
        order: BondOrder::Single,
        provenance: BondProvenance::ChemicalComponentDictionary,
    });
    data.bonds = bonds.finish();
    let bonded = Structure::new(data);
    let mut options = SaltBridgeOptions::new(4.0, SpatialBackend::BruteForce);
    assert_eq!(count(&bonded, options), Ok(0));
    options.include_bonded_pairs = true;
    assert_eq!(count(&bonded, options), Ok(1));
}

const CELL: &str = "_cell.length_a 10\n_cell.length_b 10\n_cell.length_c 10\n\
_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\n";

#[test]
fn periodic_salt_bridges_use_the_minimum_image_distance() {
    // Anion at x=5, cation at x=11.5 in a 10 Å cell: 6.5 Å directly, 3.5 Å
    // through the boundary (image at 1.5).
    let wrapped = pair(CELL, (1, 2), 11.5);
    let mut options = SaltBridgeOptions::new(4.0, SpatialBackend::BruteForce);
    assert_eq!(count(&wrapped, options), Ok(0));
    options.periodic = true;
    let Ok(table) = salt_bridges_with_options(&wrapped, options, &ExecutionContext::default())
    else {
        panic!("periodic search");
    };
    assert_eq!(table.len(), 1);
    let row = table.row(0).expect("one bridge");
    assert!((row.distance - 3.5).abs() < 1e-4, "got {}", row.distance);
}

#[test]
fn periodic_requests_reject_missing_and_placeholder_cells() {
    let mut options = SaltBridgeOptions::new(4.0, SpatialBackend::BruteForce);
    options.periodic = true;
    assert_eq!(
        count(&pair("", (1, 2), 6.5), options),
        Err(SaltBridgeError::MissingCell)
    );
    let placeholder = "_cell.length_a 1\n_cell.length_b 1\n_cell.length_c 1\n\
_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\n";
    let cube = pair(placeholder, (1, 2), 6.5);
    assert_eq!(count(&cube, options), Err(SaltBridgeError::PlaceholderCell));
    options.allow_placeholder_cell = true;
    assert!(count(&cube, options).is_ok());
}
