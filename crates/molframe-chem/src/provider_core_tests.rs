use super::*;
use crate::{ComponentKind, PolymerLinkPolicy, apply_component_chemistry};
use molframe_bench::{Sample, structure};
use molframe_core::diagnostic::Code;
use molframe_core::io::InputBuffer;
use std::fmt::Write as _;

fn get(provider: &dyn ComponentProvider, id: &str) -> Arc<Component> {
    provider
        .get(id)
        .expect("lookup works")
        .unwrap_or_else(|| panic!("{id} is bundled"))
}

#[test]
fn every_bundled_block_lowers_and_is_indexed() {
    let input = InputBuffer::from_bytes(include_bytes!("../data/core-components.cif").to_vec());
    let (full, _) = read_ccd(&input, DictionaryVersion::new("core")).expect("core set lowers");
    let core = CoreProvider;
    assert!(core.component_count() > 300);
    let lazy = LazyCifProvider::new(input, DictionaryVersion::new("core"));
    assert_eq!(lazy.block_count(), core.component_count());
    for id in ["ALA", "HOH", "ZN", "HEM", "NAG", "ATP", "DT", "PSU"] {
        assert_eq!(
            get(&core, id).atoms.len(),
            full.get(id).expect("works").expect("present").atoms.len()
        );
    }
}

#[test]
fn bundled_chemistry_matches_the_dictionary() {
    let core = CoreProvider;
    let phe = get(&core, "PHE");
    assert_eq!(phe.kind, ComponentKind::AminoAcid);
    let aromatic = phe.atoms.iter().filter(|atom| atom.aromatic).count();
    assert_eq!(aromatic, 6);
    assert!(phe.atoms.iter().any(|atom| atom.element.symbol() == "H"));
    assert_eq!(get(&core, "ZN").atoms[0].charge, 2);
    assert_eq!(get(&core, "ZN").kind, ComponentKind::Ion);
    assert_eq!(get(&core, "HOH").kind, ComponentKind::Solvent);
    assert_eq!(get(&core, "DOD").kind, ComponentKind::Solvent);
    assert_eq!(get(&core, "NAG").kind, ComponentKind::Saccharide);
    assert_eq!(get(&core, "DT").kind, ComponentKind::Nucleotide);
    assert_eq!(core.version().as_str(), CORE_DICTIONARY_VERSION);
    let ala = get(&core, "ALA");
    assert_eq!(
        ala.ideal_coordinates.as_ref().map(|ideal| ideal.len()),
        Some(ala.atoms.len())
    );
    assert!(get(&core, "GLY").ideal_coordinates.is_none());
    assert!(ala.model_coordinates.is_none());
}

#[test]
fn a_component_outside_the_core_set_is_absent_not_guessed() {
    assert!(CoreProvider.get("ZZZZ").expect("lookup works").is_none());
}

#[test]
fn the_core_set_covers_every_residue_of_the_bundled_proteins() {
    for sample in [Sample::Tiny, Sample::Medium] {
        let report = apply_component_chemistry(
            &structure(sample),
            &CoreProvider,
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
}

#[test]
fn layers_answer_in_priority_order_and_name_every_version() {
    let fake = |version: &str, name: &str| {
        let component = get(&CoreProvider, "ALA");
        let mut component = (*component).clone();
        component.name = name.into();
        let provider =
            MemoryProvider::new(DictionaryVersion::new(version), [component]).expect("unique");
        Arc::new(provider) as Arc<dyn ComponentProvider>
    };
    let layered = LayeredProvider::new(vec![fake("mine", "local alanine"), Arc::new(CoreProvider)]);
    assert_eq!(get(&layered, "ALA").name.as_ref(), "local alanine");
    assert_eq!(get(&layered, "GLY").name.to_lowercase(), "glycine");
    assert_eq!(
        layered.version().as_str(),
        format!("mine+{CORE_DICTIONARY_VERSION}")
    );
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

fn first_chain_kind(components: &[&str]) -> molframe_core::PolymerKind {
    let structure = nucleic_chain(components);
    let report = apply_component_chemistry(&structure, &CoreProvider, PolymerLinkPolicy::Disabled)
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
