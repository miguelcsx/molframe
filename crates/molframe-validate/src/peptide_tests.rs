use super::cis_peptides;
use molframe_chem::PolymerAtomRole;
use molframe_core::bond::{BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::index::AtomIndex;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_core::{AnnotationColumn, AtomAnnotation, BondOrder};

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => annotated(&structure),
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn annotated(structure: &Structure) -> Structure {
    let mut data = structure.data().clone();
    let roles = [
        PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
        PolymerAtomRole::PROTEIN_NITROGEN,
        PolymerAtomRole::PROTEIN_ALPHA_CARBON,
    ];
    let _ = data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_values(roles.into_iter().map(PolymerAtomRole::code).collect())
                .expect("small annotation column"),
        ),
    );
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: AtomIndex::new(1),
        atom_b: AtomIndex::new(2),
        order: BondOrder::Single,
        provenance: BondProvenance::User,
    });
    data.bonds = bonds.finish();
    Structure::new(data)
}

// CA(i) and CA(i+1) on the same side of the C-N bond give ω ≈ 0: a cis bond.
const CIS: &str = "\
ATOM 1 C CA GLY A 1 -0.5 1 0\n\
ATOM 2 C C GLY A 1 0 0 0\n\
ATOM 3 N N GLY A 2 1.33 0 0\n\
ATOM 4 C CA GLY A 2 1.83 1 0\n";

// The same, with CA(i+1) reflected to the far side: ω ≈ 180, a trans bond.
const TRANS: &str = "\
ATOM 1 C CA GLY A 1 -0.5 1 0\n\
ATOM 2 C C GLY A 1 0 0 0\n\
ATOM 3 N N GLY A 2 1.33 0 0\n\
ATOM 4 C CA GLY A 2 1.83 -1 0\n";

#[test]
fn a_cis_bond_is_flagged_against_the_following_residue() {
    let findings = cis_peptides(&structure(&format!("{HEADER}{CIS}")), 30.0)
        .unwrap_or_else(|error| panic!("cis validation failed: {error}"));
    assert_eq!(findings.len(), 1);
    assert_eq!(findings[0].residue.get(), 1);
    assert!(
        findings[0].omega.abs() < 1e-3,
        "ω was {}",
        findings[0].omega
    );
}

#[test]
fn a_trans_bond_is_not_flagged() {
    let findings = cis_peptides(&structure(&format!("{HEADER}{TRANS}")), 30.0)
        .unwrap_or_else(|error| panic!("cis validation failed: {error}"));
    assert!(findings.is_empty());
}
