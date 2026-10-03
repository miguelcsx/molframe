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
            Ss::AlphaHelix,
            Ss::AlphaHelix,
            Ss::AlphaHelix,
            Ss::Unknown,
            Ss::Strand,
            Ss::Strand,
            Ss::Strand
        ]
    );
    let sources = structure.data().secondary_source.to_vec();
    assert_eq!(sources[0], molframe_core::SecondarySource::None);
    assert_eq!(sources[1], molframe_core::SecondarySource::File);
    assert_eq!(sources[5], molframe_core::SecondarySource::File);
}

#[test]
fn helix_types_and_classes_name_the_helix() {
    use super::conf_state;
    assert_eq!(conf_state(Some("HELX_RH_AL_P"), None), Some(Ss::AlphaHelix));
    assert_eq!(
        conf_state(Some("HELX_RH_3T_P"), None),
        Some(Ss::ThreeTenHelix)
    );
    assert_eq!(conf_state(Some("HELX_RH_PI_P"), None), Some(Ss::PiHelix));
    assert_eq!(conf_state(Some("helx_p"), None), Some(Ss::AlphaHelix));
    assert_eq!(conf_state(Some("HELX_P"), Some(1)), Some(Ss::AlphaHelix));
    assert_eq!(conf_state(Some("HELX_P"), Some(3)), Some(Ss::PiHelix));
    assert_eq!(conf_state(Some("HELX_P"), Some(5)), Some(Ss::ThreeTenHelix));
    assert_eq!(conf_state(Some("HELX_P"), Some(7)), Some(Ss::OtherHelix));
    assert_eq!(conf_state(Some("HELX_LH_PP_P"), None), Some(Ss::OtherHelix));
    assert_eq!(conf_state(Some("TURN_P"), None), Some(Ss::Turn));
    assert_eq!(conf_state(Some("STRN"), None), None);
}

#[test]
fn a_three_ten_helix_row_reaches_the_structure() {
    let text = source().replace("HELX_P HELX_P1", "HELX_RH_3T_P HELX_RH_3T_P1");
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    assert_eq!(structure.data().secondary_structure[2], Ss::ThreeTenHelix);
}
