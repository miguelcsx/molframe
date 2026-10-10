//! Fixtures that build bonded aromatic structures from explicit coordinates.

use std::fmt::Write;

use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::index::AtomIndex;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

/// One residue: a name, a sequence number and its atom positions.
pub(crate) struct Residue {
    pub(crate) name: &'static str,
    pub(crate) atoms: Vec<[f64; 3]>,
}

/// Six ring positions of radius 1.4 Å about `centre`, spanned by the unit
/// vectors `u` and `v`.
pub(crate) fn hexagon(centre: [f64; 3], u: [f64; 3], v: [f64; 3]) -> Vec<[f64; 3]> {
    (0..6)
        .map(|k| {
            let angle = f64::from(k) * std::f64::consts::FRAC_PI_3;
            std::array::from_fn(|axis| {
                centre[axis] + 1.4 * (angle.cos() * u[axis] + angle.sin() * v[axis])
            })
        })
        .collect()
}

/// Builds a structure with every atom aromatic and the listed global bonds.
pub(crate) fn structure(residues: &[Residue], bonds: &[(u32, u32)]) -> Structure {
    let mut source = String::from(
        "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let mut serial = 0;
    for (number, residue) in residues.iter().enumerate() {
        for position in &residue.atoms {
            serial += 1;
            let _ = writeln!(
                source,
                "ATOM {serial} C C{serial} {} A {} {} {} {}",
                residue.name,
                number + 1,
                position[0],
                position[1],
                position[2]
            );
        }
    }
    let input = InputBuffer::from_bytes(source.into_bytes());
    let Ok((structure, _)) = molframe_cif::read(&input, &ReadOptions::new()) else {
        panic!("fixture must parse");
    };
    let mut data = structure.data().clone();
    let Ok(column) = AnnotationColumn::from_values(vec![true; serial]) else {
        panic!("small annotation column");
    };
    data.annotations.insert(
        molframe_core::AROMATIC_ATOM_ANNOTATION,
        AtomAnnotation::Boolean(column),
    );
    let mut table = BondTableBuilder::new();
    for &(a, b) in bonds {
        table.push(BondRecord {
            atom_a: AtomIndex::new(a),
            atom_b: AtomIndex::new(b),
            order: BondOrder::Single,
            provenance: BondProvenance::ChemicalComponentDictionary,
        });
    }
    data.bonds = table.finish();
    Structure::new(data)
}

/// Bonds closing a cycle over consecutive atoms starting at `first`.
pub(crate) fn cycle_bonds(first: u32, size: u32) -> Vec<(u32, u32)> {
    (0..size)
        .map(|k| (first + k, first + (k + 1) % size))
        .collect()
}
