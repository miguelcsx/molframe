use super::nonplanar_aromatic_rings;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

/// Bonds closing a cycle over six consecutive atoms starting at `first`.
fn ring_bonds(first: u32) -> Vec<(u32, u32)> {
    (0..6).map(|k| (first + k, first + (k + 1) % 6)).collect()
}

fn aromatic_structure(source: &str) -> Structure {
    bonded_structure(source, &ring_bonds(0))
}

fn bonded_structure(source: &str, bonds: &[(u32, u32)]) -> Structure {
    let structure = structure(source);
    let mut data = structure.data().clone();
    data.annotations.insert(
        molframe_core::AROMATIC_ATOM_ANNOTATION,
        molframe_core::AtomAnnotation::Boolean(
            molframe_core::AnnotationColumn::from_values(vec![
                true;
                structure.atom_count() as usize
            ])
            .expect("small annotation column"),
        ),
    );
    let mut table = molframe_core::BondTableBuilder::new();
    for &(a, b) in bonds {
        table.push(molframe_core::BondRecord {
            atom_a: molframe_core::AtomIndex::new(a),
            atom_b: molframe_core::AtomIndex::new(b),
            order: molframe_core::BondOrder::Single,
            provenance: molframe_core::BondProvenance::ChemicalComponentDictionary,
        });
    }
    data.bonds = table.finish();
    Structure::new(data)
}

fn flags(structure: &Structure) -> Vec<super::PlanarityFlag> {
    match nonplanar_aromatic_rings(
        structure,
        super::PlanarityOptions {
            maximum_deviation: 0.05,
            plane_fit: molframe_geom::EigenOptions::standard(),
        },
    ) {
        Ok(flags) => flags,
        Err(error) => panic!("planarity fixture failed: {error}"),
    }
}

// A phenylalanine ring laid flat in the z = 0 plane.
const FLAT_RING: &str = "\
ATOM 1 C CG PHE A 1 0.0 0.0 0\n\
ATOM 2 C CD1 PHE A 1 1.4 0.0 0\n\
ATOM 3 C CD2 PHE A 1 0.7 1.2 0\n\
ATOM 4 C CE1 PHE A 1 2.1 1.2 0\n\
ATOM 5 C CE2 PHE A 1 1.4 2.4 0\n\
ATOM 6 C CZ PHE A 1 0.0 1.2 0\n";

#[test]
fn a_flat_aromatic_ring_is_not_flagged() {
    let flags = flags(&aromatic_structure(&format!("{HEADER}{FLAT_RING}")));
    assert!(flags.is_empty(), "a flat ring should pass: {flags:?}");
}

#[test]
fn a_puckered_ring_is_flagged() {
    // The same ring with CZ lifted 1 Å out of the plane.
    let puckered = "\
ATOM 1 C CG PHE A 1 0.0 0.0 0\n\
ATOM 2 C CD1 PHE A 1 1.4 0.0 0\n\
ATOM 3 C CD2 PHE A 1 0.7 1.2 0\n\
ATOM 4 C CE1 PHE A 1 2.1 1.2 0\n\
ATOM 5 C CE2 PHE A 1 1.4 2.4 0\n\
ATOM 6 C CZ PHE A 1 0.0 1.2 1.0\n";
    let flags = flags(&aromatic_structure(&format!("{HEADER}{puckered}")));
    assert_eq!(flags.len(), 1);
    assert_eq!(flags[0].residue.get(), 0);
    assert!(flags[0].deviation > 0.05);
}

#[test]
fn a_non_aromatic_residue_is_ignored() {
    let source = format!(
        "{HEADER}\
ATOM 1 C CA ALA A 1 0.0 0.0 0\n\
ATOM 2 C CB ALA A 1 1.5 0.0 0\n"
    );
    let flags = flags(&structure(&source));
    assert!(flags.is_empty());
}

#[test]
fn invalid_plane_fit_controls_are_returned() {
    let error = nonplanar_aromatic_rings(
        &aromatic_structure(&format!("{HEADER}{FLAT_RING}")),
        super::PlanarityOptions {
            maximum_deviation: 0.05,
            plane_fit: molframe_geom::EigenOptions {
                relative_tolerance: 0.0,
                maximum_sweeps: 1,
            },
        },
    )
    .expect_err("invalid plane-fit controls must be returned");
    assert!(matches!(
        error,
        super::PlanarityError::Geometry(molframe_geom::EigenError::InvalidOptions)
    ));
}

/// Six ring atoms of radius 1.4 Å about `centre`, spanned by unit vectors `u`, `v`.
fn ring_atoms(
    first_serial: u32,
    residue: &str,
    centre: [f64; 3],
    spans: ([f64; 3], [f64; 3]),
    lift_last: f64,
) -> String {
    let (u, v) = spans;
    (0..6)
        .map(|k| {
            let serial = first_serial + k;
            let angle = f64::from(k) * std::f64::consts::FRAC_PI_3;
            let p: Vec<f64> = (0..3)
                .map(|axis| centre[axis] + 1.4 * (angle.cos() * u[axis] + angle.sin() * v[axis]))
                .collect();
            format!(
                "ATOM {serial} C C{serial} {residue} A 1 {} {} {}\n",
                p[0],
                p[1],
                p[2] + if k == 5 { lift_last } else { 0.0 }
            )
        })
        .collect::<Vec<_>>()
        .concat()
}

fn biphenyl(twist_degrees: f64, lift_second_ring_atom: f64) -> Structure {
    let angle = twist_degrees.to_radians();
    let twisted = [0.0, angle.cos(), angle.sin()];
    let mut source = HEADER.to_string();
    let first = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    source += &ring_atoms(1, "BIP", [0.0; 3], first, 0.0);
    source += &ring_atoms(
        7,
        "BIP",
        [4.3, 0.0, 0.0],
        ([1.0, 0.0, 0.0], twisted),
        lift_second_ring_atom,
    );
    let mut bonds = ring_bonds(0);
    bonds.extend(ring_bonds(6));
    bonds.push((0, 6));
    bonded_structure(&source, &bonds)
}

#[test]
fn biphenyl_with_rings_twisted_sixty_degrees_is_not_flagged() {
    let flags = flags(&biphenyl(60.0, 0.0));
    assert!(flags.is_empty(), "each ring is flat on its own: {flags:?}");
}

#[test]
fn only_the_puckered_ring_of_a_biphenyl_is_flagged() {
    let flags = flags(&biphenyl(60.0, 1.0));
    assert_eq!(flags.len(), 1, "{flags:?}");
    assert_eq!(flags[0].first_atom.get(), 6);
}
