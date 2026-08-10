use super::completeness;
use pdbiox_core::contract::Namespace;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;

// Entity A-G-S is declared; only A (position 1) and S (position 3) are modelled,
// so G at position 2 is missing.
const SOURCE: &str = r"data_seq
_entity.id 1
_entity.type polymer
_struct_asym.id A
_struct_asym.entity_id 1
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 ALA
1 2 GLY
1 3 SER
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA ALA A X 1 0 0 0
2 C CA SER A X 3 1 0 0
";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn a_gap_in_the_middle_is_reported_against_its_position() {
    let Ok(report) = completeness(&structure(), Namespace::Label) else {
        panic!("label namespace is available");
    };
    assert_eq!(report.len(), 1);
    let chain = &report[0];
    assert_eq!(chain.observed, 2);
    assert_eq!(chain.canonical, 3);
    assert_eq!(chain.missing.len(), 1);
    assert_eq!(chain.missing[0].canonical_position, 2);
}
