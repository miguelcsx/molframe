use molframe_core::SecondaryStructure as Ss;
use molframe_core::io::{InputBuffer, ReadOptions};

const ATOM_HEAD: &str = "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn source() -> String {
    let mut text = String::from("data_ss\n");
    text.push_str(ATOM_HEAD);
    text.extend((1..=8).map(|index| {
        format!(
            "ATOM {index} C CA GLY A {index} {index} {x} 0 0\n",
            x = f64::from(index) * 3.8
        )
    }));
    text.push_str(
        "loop_\n_struct_conf.conf_type_id\n_struct_conf.id\n\
_struct_conf.beg_label_comp_id\n_struct_conf.beg_label_asym_id\n\
_struct_conf.beg_label_seq_id\n_struct_conf.end_label_comp_id\n\
_struct_conf.end_label_asym_id\n_struct_conf.end_label_seq_id\n\
HELX_P HELX_P1 GLY A 2 GLY A 4\n\
loop_\n_struct_sheet_range.sheet_id\n_struct_sheet_range.id\n\
_struct_sheet_range.beg_label_comp_id\n_struct_sheet_range.beg_label_asym_id\n\
_struct_sheet_range.beg_label_seq_id\n_struct_sheet_range.end_label_comp_id\n\
_struct_sheet_range.end_label_asym_id\n_struct_sheet_range.end_label_seq_id\n\
S1 1 GLY A 6 GLY A 8\n",
    );
    text
}

#[test]
fn helix_and_sheet_ranges_label_exactly_the_residues_they_name() {
    let input = InputBuffer::from_bytes(source().into_bytes());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let states = structure.data().secondary_structure.to_vec();
    assert_eq!(
        states,
        [
            Ss::Unknown,
            Ss::Helix,
            Ss::Helix,
            Ss::Helix,
            Ss::Unknown,
            Ss::Strand,
            Ss::Strand,
            Ss::Strand
        ]
    );
}
