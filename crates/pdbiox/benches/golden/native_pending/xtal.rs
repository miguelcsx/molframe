use std::fmt::Debug;

use criterion::{BenchmarkGroup, Throughput, black_box};
use pdbiox::{
    AssemblyExt, InputBuffer, ModelIndex, ReadOptions, SpatialBackend, Structure, SymmetryExt,
};

trait BenchRequired<T> {
    fn required(self, context: &str) -> T;
}

impl<T, E: Debug> BenchRequired<T> for Result<T, E> {
    fn required(self, context: &str) -> T {
        match self {
            Ok(value) => value,
            Err(error) => panic!("{context}: {error:?}"),
        }
    }
}

impl<T> BenchRequired<T> for Option<T> {
    fn required(self, context: &str) -> T {
        match self {
            Some(value) => value,
            None => panic!("{context}"),
        }
    }
}

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

pub(super) fn register(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    bench_gw_017(group);
    bench_gw_018(group);
    bench_gw_019(group);
    bench_gw_020(group);
    bench_gw_021(group);
}

fn bench_gw_017(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = assembly_structure();
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-017", |b| {
        b.iter(|| {
            let assemblies = structure.assembly_set().required("GW-017 metadata missing");
            let view = structure.assembly("1").required("GW-017 view failed");
            black_box((
                assemblies.assemblies().count(),
                view.instance_count(),
                view.atoms().count(),
            ));
        });
    });
}

fn bench_gw_018(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = assembly_structure();
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-018", |b| {
        b.iter(|| {
            let view = structure.assembly("1").required("GW-018 view failed");
            let materialized = view.materialize().required("GW-018 materialisation failed");
            let contacts =
                pdbiox::analysis::atom_contacts(&materialized, 3.0, SpatialBackend::BruteForce)
                    .required("GW-018 contacts failed");
            black_box((materialized.chain_count(), contacts.len()));
        });
    });
}

fn bench_gw_019(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = assembly_structure();
    group.bench_function("GW-019", |b| {
        b.iter(|| {
            let view = structure.assembly("1").required("GW-019 view failed");
            let lazy = view
                .neighbors(ModelIndex::new(0), 3.0, SpatialBackend::BruteForce)
                .required("GW-019 lazy query failed");
            let materialized = view.materialize().required("GW-019 materialisation failed");
            let eager =
                pdbiox::analysis::atom_contacts(&materialized, 3.0, SpatialBackend::BruteForce)
                    .required("GW-019 eager query failed");
            black_box((lazy.len(), eager.len()));
        });
    });
}

fn bench_gw_020(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = crystal_structure();
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-020", |b| {
        b.iter(|| {
            let neighbors = structure.crystal_neighbors(10.1).required("GW-020 failed");
            black_box(neighbors.len());
        });
    });
}

fn bench_gw_021(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = assembly_structure();
    group.bench_function("GW-021", |b| {
        b.iter(|| {
            let view = structure.assembly("1").required("GW-021 view failed");
            let materialized = view.materialize().required("GW-021 materialisation failed");
            black_box((structure.atom_count(), materialized.atom_count()));
        });
    });
}

fn assembly_structure() -> Structure {
    let input = InputBuffer::from_bytes(ASSEMBLY_CIF.as_bytes().to_vec());
    let document = pdbiox::cif::parse(&input)
        .required("assembly document parse failed")
        .0;
    let assemblies = pdbiox::lower_assemblies(&document).required("assembly lowering failed");
    let structure = pdbiox::read_bytes(
        ASSEMBLY_CIF.as_bytes().to_vec(),
        Some("assembly.cif"),
        &ReadOptions::new(),
    )
    .required("assembly structure read failed")
    .0;
    structure.with_extension(pdbiox::ASSEMBLIES_EXTENSION, assemblies)
}

fn crystal_structure() -> Structure {
    let input = InputBuffer::from_bytes(CRYSTAL_CIF.as_bytes().to_vec());
    let document = pdbiox::cif::parse(&input)
        .required("crystal document parse failed")
        .0;
    let symmetry = pdbiox::lower_symmetry(&document).required("symmetry lowering failed");
    let structure = pdbiox::read_bytes(
        CRYSTAL_CIF.as_bytes().to_vec(),
        Some("crystal.cif"),
        &ReadOptions::new(),
    )
    .required("crystal structure read failed")
    .0;
    structure.with_extension(pdbiox::SYMMETRY_EXTENSION, symmetry)
}
