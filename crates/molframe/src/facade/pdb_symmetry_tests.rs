//! A PDB `CRYST1` record carries a space-group name; reading it resolves the
//! same symmetry operations a CIF read attaches.

use crate::{ReadOptions, read_bytes, write_pdb};

fn pdb(symbol: &str) -> String {
    format!(
        "CRYST1   10.000   11.000   12.000  90.00  90.00  90.00 {symbol:<11}   4\n\
ATOM      1  N   GLY A   1       1.000   1.000   1.000  1.00 10.00           N\n\
END\n"
    )
}

fn read(text: &str) -> crate::Structure {
    match read_bytes(text.as_bytes().to_vec(), Some("t.pdb"), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("the fixture reads: {findings:?}"),
    }
}

#[cfg(feature = "crystal")]
#[test]
fn a_known_space_group_symbol_attaches_its_operations() {
    use crate::crystal::SymmetryExt;

    let structure = read(&pdb("P 21 21 21"));
    let Some(symmetry) = structure.symmetry_set() else {
        panic!("the symbol resolves to a symmetry set")
    };
    assert_eq!(symmetry.operations().len(), 4);
    assert_eq!(symmetry.international_number, Some(19));
    assert_eq!(
        structure.metadata().space_group.as_deref(),
        Some("P 21 21 21")
    );
}

#[cfg(feature = "crystal")]
#[test]
fn an_unknown_symbol_is_kept_as_text_and_implies_no_symmetry() {
    use crate::crystal::SymmetryExt;

    let structure = read(&pdb("Q 9 9 9"));
    assert!(structure.symmetry_set().is_none());
    assert_eq!(structure.metadata().space_group.as_deref(), Some("Q 9 9 9"));
}

#[test]
fn the_space_group_survives_a_pdb_round_trip() {
    let structure = read(&pdb("P 21 21 21"));
    let text = match write_pdb(&structure, &crate::formats::pdb::PdbOptions::new()) {
        Ok(text) => text,
        Err(findings) => panic!("the structure writes: {findings:?}"),
    };
    assert!(text.contains("P 21 21 21"), "{text}");
    let again = read(&text);
    assert_eq!(again.metadata().space_group.as_deref(), Some("P 21 21 21"));
}

#[cfg(feature = "crystal")]
#[test]
fn remark_350_becomes_an_assembly_the_materialiser_can_expand() {
    use crate::crystal::AssemblyExt;

    let text = "\
REMARK 350 BIOMOLECULE: 1                                             \n\
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A                           \n\
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000  \n\
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000  \n\
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000  \n\
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000       10.00000  \n\
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000  \n\
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000  \n\
ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00 10.00           N\n\
END\n";
    let structure = read(text);
    let materialised = structure
        .assembly("1")
        .expect("the biomolecule is an assembly")
        .materialize()
        .expect("it materialises");
    let materialised: crate::Structure = materialised.into();
    assert_eq!(materialised.atom_count(), 2);
    let positions: Vec<[f32; 3]> = materialised.coordinates().to_vec();
    assert_eq!(positions, [[1.0, 2.0, 3.0], [9.0, -2.0, 3.0]]);
}
