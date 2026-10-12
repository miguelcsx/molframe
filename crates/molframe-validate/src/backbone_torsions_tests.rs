use super::{RamachandranError, backbone_torsions};
use molframe_chem::PolymerAtomRole;
use molframe_core::BondOrder;
use molframe_core::bond::{BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::index::AtomIndex;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_core::{AnnotationColumn, AtomAnnotation, OptionalSymbol};

/// Five alanines built from internal coordinates with phi = -57, psi = -47 and
/// omega = 180 degrees (an ideal right-handed alpha-helix): the expected angles are
/// known from the construction, not measured by the code under test.
const HELIX: &str = "data_h\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N ALA A 1 0.000 0.000 0.000\n\
ATOM 2 C CA ALA A 1 1.458 0.000 0.000\n\
ATOM 3 C C ALA A 1 2.009 1.422 0.000\n\
ATOM 4 N N ALA A 2 1.463 2.263 -0.872\n\
ATOM 5 C CA ALA A 2 1.899 3.650 -0.974\n\
ATOM 6 C C ALA A 2 1.768 4.370 0.364\n\
ATOM 7 N N ALA A 3 0.618 4.205 1.008\n\
ATOM 8 C CA ALA A 3 0.364 4.838 2.297\n\
ATOM 9 C C ALA A 3 1.421 4.443 3.323\n\
ATOM 10 N N ALA A 4 1.711 3.149 3.398\n\
ATOM 11 C CA ALA A 4 2.704 2.639 4.337\n\
ATOM 12 C C ALA A 4 4.057 3.309 4.126\n\
ATOM 13 N N ALA A 5 4.484 3.388 2.870\n\
ATOM 14 C CA ALA A 5 5.761 4.005 2.531\n\
ATOM 15 C C ALA A 5 5.830 5.442 3.035\n";

fn helix(with_roles: bool) -> Structure {
    let input = InputBuffer::from_bytes(HELIX.as_bytes().to_vec());
    let (structure, _): (Structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    if !with_roles {
        return structure;
    }
    let mut data = structure.data().clone();
    let mut bonds = BondTableBuilder::new();
    for (carbon, nitrogen) in [(2, 3), (5, 6), (8, 9), (11, 12)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(carbon),
            atom_b: AtomIndex::new(nitrogen),
            order: BondOrder::Single,
            provenance: BondProvenance::User,
            kind: OptionalSymbol::NONE,
        });
    }
    data.bonds = bonds.finish();
    let roles = [
        PolymerAtomRole::PROTEIN_NITROGEN,
        PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
    ];
    data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_values(
                (0..15)
                    .map(|index| roles[index % roles.len()].code())
                    .collect(),
            )
            .expect("small annotation column"),
        ),
    );
    Structure::new(data)
}

#[test]
fn an_ideal_helix_measures_its_constructed_torsions_for_the_inner_residues_only() {
    let torsions =
        backbone_torsions(&helix(true)).unwrap_or_else(|error| panic!("torsions failed: {error}"));
    // the first residue has no phi and the last no psi
    assert_eq!(torsions.len(), 3, "{torsions:?}");
    for t in &torsions {
        assert!((t.phi + 57.0).abs() < 1.5, "{t:?}");
        assert!((t.psi + 47.0).abs() < 1.5, "{t:?}");
    }
    assert!(
        torsions
            .windows(2)
            .all(|w| w[0].residue.get() < w[1].residue.get())
    );
}

#[test]
fn without_polymer_roles_it_refuses_instead_of_guessing_from_atom_names() {
    assert!(matches!(
        backbone_torsions(&helix(false)),
        Err(RamachandranError::Roles(_))
    ));
}
