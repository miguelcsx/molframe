use super::polymer_kind;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::topology::PolymerKind;

/// One protein chain, one DNA chain, one water chain, one `other` polymer.
const SOURCE: &str = "data_POLY
_entry.id POLY
loop_
_entity.id
_entity.type
1 polymer
2 polymer
3 water
4 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
2 polydeoxyribonucleotide
4 other
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 1 0 0 0
ATOM 2 P P DA B 2 1 1 0 0
HETATM 3 O O HOH C 3 . 2 0 0
ATOM 4 C C UNK D 4 1 3 0 0
";

fn read(source: &str) -> molframe_core::Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    structure
}

fn kinds(structure: &molframe_core::Structure) -> Vec<Option<PolymerKind>> {
    let chains = &structure.data().topology.chains;
    chains
        .iter()
        .map(|chain| chains.polymer_kind(chain))
        .collect()
}

#[test]
fn declared_polymer_types_set_each_chains_kind() {
    assert_eq!(
        kinds(&read(SOURCE)),
        [
            Some(PolymerKind::Protein),
            Some(PolymerKind::Dna),
            Some(PolymerKind::None),
            Some(PolymerKind::Other),
        ]
    );
}

#[test]
fn a_file_without_entity_poly_keeps_an_unspecified_polymer() {
    let start = SOURCE
        .find("loop_\n_entity_poly.")
        .expect("fixture declares types");
    let end = SOURCE
        .find("loop_\n_atom_site.")
        .expect("fixture has atoms");
    let undeclared = format!("{}{}", &SOURCE[..start], &SOURCE[end..]);
    assert_eq!(kinds(&read(&undeclared))[0], Some(PolymerKind::Other));
}

#[test]
fn canonical_cif_round_trip_keeps_the_declared_kinds() {
    let structure = read(SOURCE);
    let written = crate::write_canonical(&structure).expect("canonical CIF writes");
    assert!(written.contains("_entity_poly.type"), "{written}");
    assert_eq!(kinds(&read(&written)), kinds(&structure));
}

#[test]
fn every_dictionary_type_maps_and_unmodelled_types_do_not() {
    assert_eq!(polymer_kind("polypeptide(D)"), Some(PolymerKind::Protein));
    assert_eq!(polymer_kind(" POLYPEPTIDE(L) "), Some(PolymerKind::Protein));
    assert_eq!(polymer_kind("polyribonucleotide"), Some(PolymerKind::Rna));
    assert_eq!(
        polymer_kind("polydeoxyribonucleotide/polyribonucleotide hybrid"),
        Some(PolymerKind::NucleicHybrid)
    );
    assert_eq!(
        polymer_kind("polysaccharide(L)"),
        Some(PolymerKind::Saccharide)
    );
    assert_eq!(polymer_kind("peptide nucleic acid"), None);
    assert_eq!(polymer_kind("other"), None);
}
