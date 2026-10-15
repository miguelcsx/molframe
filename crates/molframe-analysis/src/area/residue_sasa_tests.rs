use super::{ResidueSasaError, residue_areas};
use crate::cif_fixture::read_cif;
use molframe_chem::HeteroPolicy;
use molframe_chem::TableRequest;
use molframe_core::ExecutionContext;
use molframe_core::structure::Structure;
use molframe_surface::SasaOptions;

fn bondi() -> SasaOptions {
    match SasaOptions::resolve(&TableRequest::named("bondi"), molframe_bench::table_store()) {
        Ok(options) => options,
        Err(error) => panic!("bondi radii: {error}"),
    }
}
use std::fmt::Write;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

/// One atom row: (group, element, atom name, residue name, sequence number, x).
fn structure(rows: &[(&str, &str, &str, &str, u32, f32)]) -> Structure {
    let mut text = HEADER.to_owned();
    for (i, (group, element, atom, residue, seq, x)) in rows.iter().enumerate() {
        let _ = writeln!(
            text,
            "{group} {} {element} {atom} {residue} A {seq} {x} 0 0",
            i + 1
        );
    }
    read_cif(&text)
}

fn areas(structure: &Structure) -> Result<super::ResidueAreas, ResidueSasaError> {
    residue_areas(structure, &bondi(), &ExecutionContext::default())
}

#[test]
fn a_lone_atom_shows_the_whole_expanded_sphere() {
    let found = areas(&structure(&[("ATOM", "N", "N", "GLY", 1, 0.0)])).unwrap();
    let expanded = f64::from(1.55_f32 + 1.4);
    let sphere = 4.0 * std::f64::consts::PI * expanded * expanded;
    assert!((found.residues[0].area - sphere).abs() < 1e-3);
    assert_eq!(found.residues[0].atoms, 1);
}

#[test]
fn a_neighbouring_residue_buries_surface() {
    let alone = areas(&structure(&[("ATOM", "C", "C", "ALA", 1, 0.0)])).unwrap();
    let pair = areas(&structure(&[
        ("ATOM", "C", "C", "ALA", 1, 0.0),
        ("ATOM", "C", "C", "ALA", 2, 3.0),
    ]))
    .unwrap();
    assert!(pair.residues[0].area < alone.residues[0].area - 1.0);
}

#[test]
fn waters_and_ligands_are_left_out_and_shield_nothing() {
    let alone = areas(&structure(&[("ATOM", "C", "C", "ALA", 1, 0.0)])).unwrap();
    let wet = areas(&structure(&[
        ("ATOM", "C", "C", "ALA", 1, 0.0),
        ("HETATM", "O", "O", "HOH", 2, 2.5),
        ("HETATM", "C", "C1", "LIG", 3, -3.0),
    ]))
    .unwrap();
    assert_eq!(wet.residues.len(), 1);
    assert_eq!(wet.excluded_hetero_residues, 2);
    assert!((wet.residues[0].area - alone.residues[0].area).abs() < 1e-9);
}

#[test]
fn a_ligand_the_options_include_shields_its_neighbour() {
    let rows = [
        ("ATOM", "C", "C", "ALA", 1, 0.0),
        ("HETATM", "C", "C1", "LIG", 2, 3.0),
    ];
    let options = SasaOptions {
        hetero: HeteroPolicy::SOLUTE,
        ..bondi()
    };
    let context = ExecutionContext::default();
    let with = residue_areas(&structure(&rows), &options, &context).unwrap();
    let without = areas(&structure(&rows)).unwrap();
    assert_eq!(with.residues.len(), 2);
    assert_eq!(with.excluded_hetero_residues, 0);
    assert!(with.residues[0].area < without.residues[0].area - 1.0);
}

#[test]
fn a_structure_with_only_water_has_nothing_to_measure() {
    let found = areas(&structure(&[("HETATM", "O", "O", "HOH", 1, 0.0)]));
    assert!(matches!(found, Err(ResidueSasaError::NoAtoms)));
}

#[test]
fn a_residues_area_is_the_sum_of_its_atoms_areas() {
    let rows = [
        ("ATOM", "N", "N", "ALA", 1, 0.0),
        ("ATOM", "C", "CA", "ALA", 1, 1.5),
        ("ATOM", "H", "H", "ALA", 1, 0.5),
        ("ATOM", "C", "C", "GLY", 2, 3.0),
        ("ATOM", "O", "O", "GLY", 2, 4.2),
        ("HETATM", "O", "O", "HOH", 3, 2.0),
    ];
    let structure = structure(&rows);
    let options = bondi();
    let context = ExecutionContext::default();
    let found = residue_areas(&structure, &options, &context).unwrap();
    let radii = crate::sasa_atoms::measured_radii(&structure, &options).radii;
    assert!(radii[2].is_nan() && radii[5].is_nan());
    let atoms =
        molframe_surface::measure_atoms(structure.positions(), &radii, &options, &context).unwrap();
    let sums = [atoms[0] + atoms[1], atoms[3] + atoms[4]];
    assert_eq!(found.residues.len(), 2);
    for (residue, sum) in found.residues.iter().zip(sums) {
        assert_eq!(residue.area.to_bits(), sum.to_bits());
    }
    assert_eq!(found.residues[0].atoms, 2);
    assert_eq!(found.missing_radius, 0);
}
