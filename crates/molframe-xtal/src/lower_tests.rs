use super::lower_assemblies;
use molframe_cif::parse;
use molframe_core::{Code, InputBuffer};

const ASSEMBLY: &str = r"data_demo
loop_
_pdbx_struct_oper_list.id
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
1 1 0 0 0 0 1 0 0 0 0 1 0
2 1 0 0 10 0 1 0 20 0 0 1 30
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.method_details
_pdbx_struct_assembly.oligomeric_count
1 'author_defined_assembly' PISA 4
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1,2)(1-2)' 'A,B'
";

#[test]
fn all_three_categories_lower_without_losing_identifiers() {
    let document = parsed(ASSEMBLY);
    let set = match lower_assemblies(&document) {
        Ok(set) => set,
        Err(findings) => panic!("lowering failed: {findings:?}"),
    };

    assert_eq!(set.len(), 1);
    assert_eq!(set.operators().count(), 2);
    let Some(definition) = set.get("1") else {
        panic!("assembly absent")
    };
    assert_eq!(
        definition.details.as_deref(),
        Some("author_defined_assembly")
    );
    assert_eq!(definition.method.as_deref(), Some("PISA"));
    assert_eq!(definition.oligomeric, Some(4));
    assert_eq!(definition.generators.len(), 1);
    assert_eq!(
        definition.generators[0].oper_expression.combination_count(),
        4
    );
    assert_eq!(
        definition.generators[0]
            .asym_ids
            .iter()
            .map(Box::as_ref)
            .collect::<Vec<_>>(),
        ["A", "B"]
    );
    assert_eq!(
        set.operator("2")
            .map(|operator| operator.transform.translation),
        Some([10.0, 20.0, 30.0])
    );
}

#[test]
fn a_generator_cannot_name_an_absent_assembly() {
    let document = parsed(
        "data_demo\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n9 1 A\n",
    );
    let Err(findings) = lower_assemblies(&document) else {
        panic!("unknown assembly accepted")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E6013));
}

#[test]
fn a_non_orthogonal_operator_is_rejected() {
    let malformed = ASSEMBLY.replacen("1 1 0 0 0 0 1", "1 2 0 0 0 0 1", 1);
    let document = parsed(&malformed);
    let Err(findings) = lower_assemblies(&document) else {
        panic!("non-rigid operator accepted")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E6012));
}

fn parsed(text: &str) -> molframe_cif::Document {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match parse(&input) {
        Ok((document, findings)) if findings.is_empty() => document,
        Ok((_, findings)) => panic!("fixture diagnostics: {findings:?}"),
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    }
}
