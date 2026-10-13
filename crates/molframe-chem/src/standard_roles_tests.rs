use super::*;
use crate::assign_secondary_structure;
use molframe_bench::{Sample, structure};
use molframe_core::annotation::AtomAnnotation;
use molframe_core::column::Presence;

fn role_of(structure: &Structure, atom: u32) -> (PolymerAtomRole, Presence) {
    let Some(AtomAnnotation::Integer(column)) = structure
        .annotations()
        .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION)
    else {
        panic!("role column missing");
    };
    let (code, presence) = column.get(atom).expect("row exists");
    (
        PolymerAtomRole::from_code(code).expect("known bits"),
        presence,
    )
}

#[test]
fn backbone_and_beta_atoms_of_amino_acids_get_their_roles() {
    let source = structure(Sample::Tiny);
    let annotated = annotate_standard_protein_roles(&source).expect("small column");
    let mut backbone = 0;
    for atom in annotated.data().atoms() {
        let (role, presence) = role_of(&annotated, atom.index().get());
        match atom.name() {
            Some("N" | "CA" | "C" | "O") => {
                assert!(role.intersects(PolymerAtomRole::PROTEIN_BACKBONE));
                assert_eq!(presence, Presence::Present);
                backbone += 1;
            }
            Some("CB") => assert!(role.intersects(PolymerAtomRole::PROTEIN_BETA_CARBON)),
            _ => assert_eq!(role, PolymerAtomRole::UNKNOWN),
        }
    }
    assert_eq!(backbone, 4 * source.residue_count());
}

#[test]
fn ligands_and_waters_keep_unknown_roles_even_when_they_own_a_ca_or_c() {
    let source = structure(Sample::Medium);
    let annotated = annotate_standard_protein_roles(&source).expect("small column");
    for atom in annotated.data().atoms() {
        let name = atom.residue().and_then(|residue| residue.name());
        if matches!(name, Some("HEM" | "HOH")) {
            let (role, presence) = role_of(&annotated, atom.index().get());
            assert_eq!(role, PolymerAtomRole::UNKNOWN, "{name:?}");
            assert_eq!(presence, Presence::Unknown);
        }
    }
}

#[test]
fn an_existing_role_column_is_never_replaced() {
    let source = structure(Sample::Tiny);
    let first = annotate_standard_protein_roles(&source).expect("small column");
    let again = annotate_standard_protein_roles(&first).expect("idempotent");
    for atom in first.data().atoms() {
        let index = atom.index().get();
        assert_eq!(role_of(&first, index), role_of(&again, index));
    }
}

fn letters(source: &Structure) -> String {
    use molframe_core::SecondaryStructure as Ss;
    assign_secondary_structure(source)
        .iter()
        .map(|a| match a.state {
            Ss::AlphaHelix => 'H',
            Ss::ThreeTenHelix => 'G',
            Ss::PiHelix => 'I',
            Ss::PolyProline => 'P',
            Ss::Strand => 'E',
            Ss::BetaBridge => 'B',
            Ss::Turn => 'T',
            Ss::Bend => 'S',
            Ss::Coil => '-',
            _ => '?',
        })
        .collect()
}

/// One-based inclusive residue range of `letters` must hold only `expected`.
fn holds(letters: &str, first: usize, last: usize, expected: char) -> bool {
    letters[first - 1..last].chars().all(|c| c == expected)
}

#[test]
fn ubiquitin_matches_its_published_secondary_structure() {
    let line = letters(&structure(Sample::Small));
    assert!(holds(&line, 2, 7, 'E'), "beta1 {line}");
    assert!(holds(&line, 12, 16, 'E'), "beta2 {line}");
    assert!(holds(&line, 23, 34, 'H'), "alpha1 {line}");
    assert!(holds(&line, 38, 40, 'G'), "3-10 {line}");
    assert!(holds(&line, 41, 45, 'E'), "beta3 {line}");
    assert!(holds(&line, 57, 59, 'G'), "3-10 {line}");
    assert!(holds(&line, 66, 71, 'E'), "beta5 {line}");
    assert!(line[..76].chars().all(|c| c != '?'), "protein is evaluable");
    assert!(line[76..].chars().all(|c| c == '?'), "waters stay unknown");
}

#[test]
fn crambin_matches_its_published_secondary_structure() {
    let line = letters(&structure(Sample::Tiny));
    assert!(holds(&line, 7, 17, 'H'), "helix 1 {line}");
    assert!(holds(&line, 23, 29, 'H'), "helix 2 {line}");
    assert!(holds(&line, 2, 3, 'E'), "strand 1 {line}");
    assert!(holds(&line, 33, 34, 'E'), "strand 2 {line}");
    assert!(holds(&line, 42, 44, 'G'), "3-10 {line}");
}
