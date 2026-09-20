use molframe::crystal::{AssemblyExt, SymmetryExt};
use molframe::spatial::SpatialBackend;
use molframe::{InputBuffer, ModelIndex, ReadOptions, Structure};

const ASSEMBLY_CIF: &str = r"data_demo
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

const CRYSTAL_CIF: &str = r"data_crystal
_cell.length_a 10
_cell.length_b 100
_cell.length_c 100
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
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
1 C CA GLY A 1 1 A 1 0 0
loop_
_space_group_symop.id
_space_group_symop.operation_xyz
1 'x,y,z'
";

#[test]
fn gw_017_lists_assemblies_without_materialising_instances() {
    let structure = assembly_structure();
    let assemblies = structure
        .assembly_set()
        .unwrap_or_else(|| panic!("assembly metadata missing"));
    assert_eq!(
        assemblies
            .assemblies()
            .map(|assembly| assembly.id.as_ref())
            .collect::<Vec<_>>(),
        ["1"]
    );
    let view = structure
        .assembly("1")
        .unwrap_or_else(|finding| panic!("assembly view failed: {finding}"));
    assert_eq!(view.instance_count(), 2);
    assert_eq!(view.atoms().count(), 2);
    assert!(std::ptr::eq(
        view.source().positions().as_ptr(),
        structure.coordinates().as_ptr()
    ));
}

#[test]
fn gw_018_materialises_a_biological_assembly_and_finds_its_interface() {
    let structure = assembly_structure();
    let view = structure
        .assembly("1")
        .unwrap_or_else(|finding| panic!("assembly view failed: {finding}"));
    let materialized = view
        .materialize()
        .unwrap_or_else(|findings| panic!("assembly materialisation failed: {findings:?}"));
    let contacts = molframe::analysis::atom_contacts(
        &materialized,
        3.0,
        SpatialBackend::BruteForce,
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("assembly interface failed: {error}"));
    assert_eq!(materialized.chain_count(), 2);
    assert!(contacts.iter().any(|contact| contact.distance < 3.0));
}

#[test]
fn gw_019_lazy_and_materialised_assembly_interface_queries_agree() {
    let structure = assembly_structure();
    let view = structure
        .assembly("1")
        .unwrap_or_else(|finding| panic!("assembly view failed: {finding}"));
    let lazy = view
        .neighbors(
            ModelIndex::new(0),
            3.0,
            SpatialBackend::BruteForce,
            &molframe::ExecutionContext::default(),
        )
        .unwrap_or_else(|finding| panic!("lazy assembly query failed: {finding}"));
    let materialized = view
        .materialize()
        .unwrap_or_else(|findings| panic!("assembly materialisation failed: {findings:?}"));
    let eager = molframe::analysis::atom_contacts(
        &materialized,
        3.0,
        SpatialBackend::BruteForce,
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("eager assembly query failed: {error}"));
    assert_eq!(lazy.len(), eager.len());
    assert!(lazy.iter().zip(eager.iter()).all(|(left, right)| {
        (f64::from(left.distance_squared).sqrt() - f64::from(right.distance)).abs() < 1.0e-5
    }));
}

#[test]
fn gw_020_finds_crystal_contacts_across_symmetry_images() {
    let structure = crystal_structure();
    let neighbors = structure
        .collect_crystal_neighbors(10.1, &molframe::ExecutionContext::default())
        .unwrap_or_else(|finding| panic!("crystal neighbour search failed: {finding}"));
    assert_eq!(neighbors.len(), 1);
    assert_eq!(neighbors[0].lattice, [-1, 0, 0]);
    assert!((neighbors[0].distance_squared - 100.0).abs() < 1.0e-9);
}

#[test]
fn gw_021_keeps_asymmetric_unit_and_biological_assembly_conclusions_distinct() {
    let structure = assembly_structure();
    let view = structure
        .assembly("1")
        .unwrap_or_else(|finding| panic!("assembly view failed: {finding}"));
    let materialized = view
        .materialize()
        .unwrap_or_else(|findings| panic!("assembly materialisation failed: {findings:?}"));
    assert_eq!(structure.atom_count(), 2);
    assert_eq!(materialized.atom_count(), 2);
    assert_ne!(
        structure.coordinates().as_ptr(),
        materialized.positions().as_ptr()
    );
    assert_ne!(
        materialized.positions()[0][0].to_bits(),
        structure.coordinates()[0][0].to_bits()
    );
}

fn assembly_structure() -> Structure {
    let input = InputBuffer::from_bytes(ASSEMBLY_CIF.as_bytes().to_vec());
    let document = molframe::formats::cif::parse(&input)
        .unwrap_or_else(|findings| panic!("assembly document parse failed: {findings:?}"))
        .0;
    let assemblies = molframe::crystal::lower_assemblies(&document)
        .unwrap_or_else(|findings| panic!("assembly lowering failed: {findings:?}"));
    let structure = molframe::read_bytes(
        ASSEMBLY_CIF.as_bytes().to_vec(),
        Some("assembly.cif"),
        &ReadOptions::new(),
    )
    .unwrap_or_else(|findings| panic!("assembly structure read failed: {findings:?}"))
    .0;
    structure
        .engine()
        .with_extension(molframe::crystal::ASSEMBLIES_EXTENSION, assemblies)
        .into()
}

fn crystal_structure() -> Structure {
    let input = InputBuffer::from_bytes(CRYSTAL_CIF.as_bytes().to_vec());
    let document = molframe::formats::cif::parse(&input)
        .unwrap_or_else(|findings| panic!("crystal document parse failed: {findings:?}"))
        .0;
    let symmetry = molframe::crystal::lower_symmetry(&document)
        .unwrap_or_else(|findings| panic!("symmetry lowering failed: {findings:?}"));
    let structure = molframe::read_bytes(
        CRYSTAL_CIF.as_bytes().to_vec(),
        Some("crystal.cif"),
        &ReadOptions::new(),
    )
    .unwrap_or_else(|findings| panic!("crystal structure read failed: {findings:?}"))
    .0;
    structure
        .engine()
        .with_extension(molframe::crystal::SYMMETRY_EXTENSION, symmetry)
        .into()
}
