use super::{AssemblyExt, AssemblyView};
use crate::{ASSEMBLIES_EXTENSION, lower_assemblies};
use molframe_cif::{parse, read};
use molframe_core::{Code, InputBuffer, ModelIndex, ReadOptions, Structure};

pub(crate) const ENTRY: &str = r"data_demo
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 1 X 1 0 0
2 C CA GLY B 1 1 Y 0 2 0
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
T 1 0 0 10 0 1 0 0 0 0 1 0
R 0 -1 0 0 1 0 0 0 0 0 1 0
_pdbx_struct_assembly.id 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(T)(R)' 'A,B'
";

#[test]
fn assembly_instances_are_lazy_stable_and_composed_right_to_left() {
    let structure = attached(ENTRY);
    let view = match structure.assembly("1") {
        Ok(view) => view,
        Err(finding) => panic!("assembly failed: {finding}"),
    };

    assert_eq!(view.instance_count(), 2);
    assert_eq!(view.atoms().count(), 2);
    let positions = view
        .positions(ModelIndex::new(0))
        .filter_map(|(_, position)| position)
        .collect::<Vec<_>>();
    assert_near(positions[0], [10.0, 1.0, 0.0]);
    assert_near(positions[1], [8.0, 0.0, 0.0]);
    assert_eq!(
        view.chains()
            .map(|instance| instance.instance_id.get())
            .collect::<Vec<_>>(),
        [0, 1]
    );
    assert!(std::ptr::eq(
        view.source().positions().as_ptr(),
        structure.positions().as_ptr()
    ));
}

#[test]
fn generators_resolve_only_the_label_namespace() {
    let structure = attached(&ENTRY.replace("'A,B'", "'X,Y'"));
    let finding = structure.assembly("1").err();
    assert_eq!(finding.map(|finding| finding.code()), Some(Code::E6013));
}

#[test]
fn total_chain_instance_expansion_is_bounded() {
    let structure = attached(ENTRY);
    let Some(assemblies) = structure.assembly_set() else {
        panic!("assembly set absent")
    };
    let finding = AssemblyView::with_limit(&structure, assemblies, "1", 1).err();
    assert_eq!(finding.map(|finding| finding.code()), Some(Code::E6011));
}

pub(crate) fn attached(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    let assemblies = match lower_assemblies(&document) {
        Ok(assemblies) => assemblies,
        Err(findings) => panic!("assembly lowering failed: {findings:?}"),
    };
    let structure = match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("structure lowering failed: {findings:?}"),
    };
    structure.with_extension(ASSEMBLIES_EXTENSION, assemblies)
}

fn assert_near(actual: [f32; 3], expected: [f32; 3]) {
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(actual, expected)| (*actual - expected).abs() < 1e-5),
        "{actual:?} != {expected:?}"
    );
}
