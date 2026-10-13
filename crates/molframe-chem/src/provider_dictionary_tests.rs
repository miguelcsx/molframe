use super::*;
use crate::{ComponentKind, PolymerLinkPolicy, TableError, TableStore, apply_component_chemistry};
use molframe_bench::{Sample, structure};
use molframe_core::diagnostic::Code;
use molframe_core::io::InputBuffer;
use std::fmt::Write as _;

const AMINO_ACIDS: &[u8] = include_bytes!("../data/CCD-amino-acids.cif");
const SACCHARIDES: &[u8] = include_bytes!("../data/CCD-saccharides.cif");

/// A store over a fresh directory holding `entries` and a profile selecting `chosen`.
pub(super) fn fixture_store(
    label: &str,
    entries: &[(&str, &[u8])],
    chosen: Option<&str>,
) -> TableStore {
    let dir = std::env::temp_dir().join(format!("molframe-dict-{label}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::create_dir_all(&dir);
    for (name, bytes) in entries {
        let _ = std::fs::write(dir.join(format!("{name}.cif")), bytes);
    }
    if let Some(name) = chosen {
        let _ = std::fs::write(dir.join("profile.txt"), format!("components: {name}\n"));
    }
    TableStore::new(dir)
}

fn selected(store: &TableStore) -> Arc<dyn ComponentProvider> {
    profile_dictionary(store)
        .expect("selection works")
        .expect("the profile names a dictionary")
}

fn get(provider: &dyn ComponentProvider, id: &str) -> Arc<Component> {
    provider
        .get(id)
        .expect("lookup works")
        .unwrap_or_else(|| panic!("{id} is defined"))
}

fn amino_store(label: &str) -> TableStore {
    fixture_store(label, &[("fixture-ccd", AMINO_ACIDS)], Some("fixture-ccd"))
}

#[test]
fn the_profile_selects_a_dictionary_and_its_stem_is_the_version() {
    let provider = selected(&amino_store("select"));
    assert_eq!(provider.version().as_str(), "fixture-ccd");
    assert_eq!(get(provider.as_ref(), "ALA").atoms.len(), 13);
}

#[test]
fn a_profile_without_a_components_line_selects_nothing() {
    let store = fixture_store("noline", &[("fixture-ccd", AMINO_ACIDS)], None);
    assert!(profile_dictionary(&store).expect("not an error").is_none());
}

#[test]
fn a_named_entry_the_cache_lacks_is_a_typed_missing() {
    let store = fixture_store("lacking", &[], Some("absent"));
    let Err(error) = profile_dictionary(&store) else {
        panic!("the entry is missing");
    };
    assert!(matches!(
        error,
        TableError::Missing {
            kind: "components",
            ..
        }
    ));
}

#[test]
fn dictionary_chemistry_is_read_from_the_file() {
    let provider = selected(&amino_store("chem"));
    let phe = get(provider.as_ref(), "PHE");
    assert_eq!(phe.kind, ComponentKind::AminoAcid);
    assert_eq!(phe.atoms.iter().filter(|atom| atom.aromatic).count(), 6);
    assert!(phe.atoms.iter().any(|atom| atom.element.symbol() == "H"));
    assert_eq!(
        get(provider.as_ref(), "LYS")
            .atoms
            .iter()
            .map(|a| i32::from(a.charge))
            .sum::<i32>(),
        1
    );
}

#[test]
fn saccharides_are_typed_by_the_dictionary() {
    let store = fixture_store("sugar", &[("sugars", SACCHARIDES)], Some("sugars"));
    let provider = selected(&store);
    assert_eq!(
        get(provider.as_ref(), "NAG").kind,
        ComponentKind::Saccharide
    );
}

#[test]
fn a_component_outside_the_dictionary_is_absent_not_guessed() {
    let provider = selected(&amino_store("absent"));
    assert!(provider.get("ZZZZ").expect("lookup works").is_none());
}

#[test]
fn the_dictionary_covers_every_residue_of_a_protein() {
    let provider = selected(&amino_store("cover"));
    let report = apply_component_chemistry(
        &structure(Sample::Tiny),
        provider.as_ref(),
        PolymerLinkPolicy::Disabled,
    )
    .expect("provider works");
    let unknown: Vec<_> = report
        .findings
        .iter()
        .filter(|finding| finding.code() == Code::W3201)
        .collect();
    assert!(unknown.is_empty(), "{unknown:?}");
}

#[test]
fn layers_answer_in_priority_order_and_name_every_version() {
    let fixture = selected(&amino_store("layers"));
    let component = get(fixture.as_ref(), "ALA");
    let mut component = (*component).clone();
    component.name = "local alanine".into();
    let mine = MemoryProvider::new(DictionaryVersion::new("mine"), [component]).expect("unique");
    let layered = LayeredProvider::new(vec![Arc::new(mine), fixture]);
    assert_eq!(get(&layered, "ALA").name.as_ref(), "local alanine");
    assert_eq!(get(&layered, "GLY").name.to_lowercase(), "glycine");
    assert_eq!(layered.version().as_str(), "mine+fixture-ccd");
    assert!(layered.get("ZZZZ").expect("lookup works").is_none());
}

const IN_FILE: &str = "data_s\n\
loop_\n_chem_comp.id\n_chem_comp.name\n_chem_comp.type\n\
LIG 'a ligand' NON-POLYMER\nNOA 'no atoms given' NON-POLYMER\n\
LIG2 'another' NON-POLYMER\n\
loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
_chem_comp_atom.pdbx_leaving_atom_flag\n_chem_comp_atom.pdbx_stereo_config\n\
LIG C1 C 0 N N N\nLIG O1 O 0 N N N\nLIG2 N1 N 0 N N N\n\
loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n\
_chem_comp_bond.pdbx_stereo_config\nLIG C1 O1 DOUB N N\n";

#[test]
fn in_file_definitions_are_grouped_by_component() {
    let input = InputBuffer::from_bytes(IN_FILE.as_bytes().to_vec());
    let (document, _) = molframe_cif::parse(&input).expect("parses");
    let provider = EmbeddedProvider::from_document(&document, DictionaryVersion::new("in-file"))
        .expect("lowers");
    assert!(!provider.is_empty());
    let lig = get(&provider, "LIG");
    assert_eq!(lig.atoms.len(), 2);
    assert_eq!(lig.bonds.len(), 1);
    assert_eq!(get(&provider, "LIG2").atoms.len(), 1);
    assert!(get(&provider, "LIG2").bonds.is_empty());
}

#[test]
fn a_listed_component_with_no_atom_rows_stays_unknown() {
    let input = InputBuffer::from_bytes(IN_FILE.as_bytes().to_vec());
    let (document, _) = molframe_cif::parse(&input).expect("parses");
    let provider = EmbeddedProvider::from_document(&document, DictionaryVersion::new("in-file"))
        .expect("lowers");
    assert!(provider.get("NOA").expect("lookup works").is_none());
}

#[test]
fn a_file_without_embedded_atoms_defines_nothing() {
    let input = InputBuffer::from_bytes(
        b"data_s\nloop_\n_chem_comp.id\n_chem_comp.name\n_chem_comp.type\nLIG x NON-POLYMER\n"
            .to_vec(),
    );
    let (document, _) = molframe_cif::parse(&input).expect("parses");
    let provider =
        EmbeddedProvider::from_document(&document, DictionaryVersion::new("v")).expect("lowers");
    assert!(provider.is_empty());
}

fn nucleic_chain(components: &[&str]) -> molframe_core::Structure {
    let mut text = String::from(
        "data_na\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    for (index, component) in components.iter().enumerate() {
        let number = index + 1;
        let x = number * 6;
        writeln!(text, "ATOM {number} C C1' {component} A {number} {x} 0 0").expect("string write");
    }
    let buffer = InputBuffer::from_bytes(text.into_bytes());
    molframe_cif::read(&buffer, &molframe_core::io::ReadOptions::new())
        .expect("reads")
        .0
}

/// A dictionary of nucleotides: the 2'-oxygen is what separates RNA from DNA.
fn nucleotide_dictionary() -> String {
    let mut text = String::new();
    for (id, kind, oxygen) in [
        ("A", "RNA LINKING", true),
        ("G", "RNA LINKING", true),
        ("U", "RNA LINKING", true),
        ("C", "RNA LINKING", true),
        ("PSU", "RNA LINKING", true),
        ("DA", "DNA LINKING", false),
        ("DG", "DNA LINKING", false),
        ("DT", "DNA LINKING", false),
        ("DC", "DNA LINKING", false),
    ] {
        let _ = write!(
            text,
            "data_{id}\n_chem_comp.id {id}\n_chem_comp.name {id}\n_chem_comp.type '{kind}'\n\
loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_leaving_atom_flag\n\
_chem_comp_atom.pdbx_stereo_config\n{id} C1' C 0 N N N\n{id} C2' C 0 N N N\n"
        );
        if oxygen {
            let _ = writeln!(text, "{id} O2' O 0 N N N");
        }
    }
    text
}

fn first_chain_kind(components: &[&str]) -> molframe_core::PolymerKind {
    let structure = nucleic_chain(components);
    let dictionary = nucleotide_dictionary().into_bytes();
    let store = fixture_store(
        "nucleic",
        &[("nucleotides", &dictionary)],
        Some("nucleotides"),
    );
    let provider = selected(&store);
    let report =
        apply_component_chemistry(&structure, provider.as_ref(), PolymerLinkPolicy::Disabled)
            .expect("provider works");
    let chain = report.structure.data().chains().next().expect("one chain");
    chain.polymer_kind()
}

#[test]
fn nucleic_chains_are_told_apart_by_their_two_prime_oxygen() {
    use molframe_core::PolymerKind;
    assert_eq!(
        first_chain_kind(&["A", "G", "U", "C", "PSU"]),
        PolymerKind::Rna
    );
    assert_eq!(
        first_chain_kind(&["DA", "DG", "DT", "DC"]),
        PolymerKind::Dna
    );
    assert_eq!(
        first_chain_kind(&["DA", "DG", "A", "U"]),
        PolymerKind::NucleicHybrid
    );
}
