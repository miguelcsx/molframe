use super::{assemblies_cif, ncs_cif};
use crate::{lower_assemblies, lower_ncs};
use molframe_cif::{Document, parse};
use molframe_core::InputBuffer;

const SOURCE: &str = r"data_demo
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0
2 'crystal symmetry operation' 1 0 0 10.5 0 1 0 20 0 0 1 30.25
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.method_details
_pdbx_struct_assembly.oligomeric_count
1 'author_defined_assembly' PISA 4
2 'software_defined_assembly' ? 2
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1,2)(1-2)' 'A,B'
2 1 A
loop_
_struct_ncs_oper.id
_struct_ncs_oper.code
_struct_ncs_oper.details
_struct_ncs_oper.matrix[1][1]
_struct_ncs_oper.matrix[1][2]
_struct_ncs_oper.matrix[1][3]
_struct_ncs_oper.vector[1]
_struct_ncs_oper.matrix[2][1]
_struct_ncs_oper.matrix[2][2]
_struct_ncs_oper.matrix[2][3]
_struct_ncs_oper.vector[2]
_struct_ncs_oper.matrix[3][1]
_struct_ncs_oper.matrix[3][2]
_struct_ncs_oper.matrix[3][3]
_struct_ncs_oper.vector[3]
1 given ? 1 0 0 0 0 1 0 0 0 0 1 0
2 generate 'second copy' 0 -1 0 1.5 1 0 0 2.5 0 0 1 3.5
";

fn parsed(text: &str) -> Document {
    match parse(&InputBuffer::from_bytes(text.as_bytes().to_vec())) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    }
}

#[test]
fn assemblies_survive_a_text_round_trip() {
    let first = lower_assemblies(&parsed(SOURCE)).expect("source lowers");
    let text = assemblies_cif(&first);
    let again = lower_assemblies(&parsed(&format!("data_x\n{text}"))).expect("output lowers");

    assert_eq!(again.len(), first.len());
    assert_eq!(again.operators().count(), 2);
    for (a, b) in first.assemblies().zip(again.assemblies()) {
        assert_eq!(a, b);
    }
    for (a, b) in first.operators().zip(again.operators()) {
        assert_eq!(a, b);
    }
}

#[test]
fn ncs_operators_survive_a_text_round_trip() {
    let first = lower_ncs(&parsed(SOURCE)).expect("source lowers");
    let text = ncs_cif(&first);
    let again = lower_ncs(&parsed(&format!("data_x\n{text}"))).expect("output lowers");

    assert_eq!(again.len(), 2);
    for (a, b) in first.operators().zip(again.operators()) {
        assert_eq!(a, b);
    }
}

#[test]
fn an_empty_set_writes_nothing() {
    assert!(assemblies_cif(&crate::AssemblySet::default()).is_empty());
    assert!(ncs_cif(&crate::NcsSet::default()).is_empty());
}
