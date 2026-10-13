use super::*;
use molframe_bench::{Sample, structure};
use molframe_core::io::InputBuffer;

const DEFINED: &str = "data_s\n\
loop_\n_chem_comp.id\n_chem_comp.name\n_chem_comp.type\nLIG 'a ligand' NON-POLYMER\n\
loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
_chem_comp_atom.pdbx_leaving_atom_flag\n_chem_comp_atom.pdbx_stereo_config\n\
LIG C1 C 0 N N N\nLIG O1 O 0 N N N\n";

fn document(text: &str) -> Document {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    molframe_cif::parse(&input).expect("parses").0
}

#[test]
fn a_document_without_definitions_yields_no_provider() {
    let none = embedded_components(&document("data_s\n_entry.id X\n")).expect("lowers");
    assert!(none.is_none());
}

#[test]
fn stack_answers_in_file_then_user_then_core() {
    let embedded = embedded_components(&document(DEFINED))
        .expect("lowers")
        .expect("defines LIG");
    let plain = structure(Sample::Tiny);
    let with = plain.with_extension(EMBEDDED_COMPONENTS_EXTENSION, embedded);

    let without = resolution_stack(&plain, None);
    assert!(without.get("LIG").expect("lookup").is_none());
    assert!(without.get("ALA").expect("lookup").is_some());

    let stack = resolution_stack(&with, None);
    assert_eq!(
        stack
            .get("LIG")
            .expect("lookup")
            .expect("in file")
            .atoms
            .len(),
        2
    );
    assert!(stack.get("ALA").expect("lookup").is_some());
    assert!(stack.version().as_str().starts_with("embedded+"));
}

const ENTRY_ROWS: &str = "data_s\n\
loop_\n_chem_comp.id\n_chem_comp.name\n_chem_comp.type\n\
LIG 'a ligand' NON-POLYMER\nLYS 'lysine, uncharged here' 'L-PEPTIDE LINKING'\n\
loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_stereo_config\n\
LIG C1 C N N\nLIG O1 O N N\nLYS NZ N N N\n";

#[test]
fn definitions_without_charges_rank_below_the_dictionaries() {
    let embedded = embedded_components(&document(ENTRY_ROWS))
        .expect("entry-style rows lower")
        .expect("defines components");
    assert!(!embedded.states_charges("LIG"));
    let structure = structure(Sample::Tiny).with_extension(EMBEDDED_COMPONENTS_EXTENSION, embedded);
    let stack = resolution_stack(&structure, None);
    assert_eq!(
        stack
            .get("LIG")
            .expect("lookup")
            .expect("in file")
            .atoms
            .len(),
        2
    );
    let lysine = stack.get("LYS").expect("lookup").expect("core answers");
    assert!(
        lysine.atoms.len() > 1,
        "the core definition wins over a charge-less one"
    );
}
