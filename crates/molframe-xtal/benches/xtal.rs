//! Criterion coverage for crystal catalogue and symmetry workflows.

use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, Criterion, black_box};
use molframe_bench::{Sample, input, structure};
use molframe_core::ModelIndex;
use molframe_core::execution::ExecutionContext;
use molframe_core::structure::UnitCell;
use molframe_xtal::{
    AffineTransform, CellTransform, DensityMap, MapBoundary, MrcBlockOptions, MrcBlockReader,
    NcsSet, ReflectionColumn, ReflectionColumnType, ReflectionTable, ReflectionValue,
    collect_crystal_neighbors, lower_assemblies, lower_ncs, lower_structure_factor_cif,
    lower_symmetry, read_cube, read_dx, read_mtz, space_group_by_hall, space_group_setting,
    space_group_settings, write_mtz, write_structure_factor_cif,
};
use std::fmt::{Debug, Write};
use std::io::Cursor;

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

fn bench_space_groups(c: &mut Criterion) {
    let mut group = c.benchmark_group("xtal_space_groups");
    group.bench_function("hall_lookup", |b| {
        b.iter(|| black_box(space_group_by_hall(" P 1 ")));
    });
    group.bench_function("number_lookup", |b| {
        b.iter(|| black_box(space_group_setting(1)));
    });
    group.bench_function("settings_for_type", |b| {
        b.iter(|| black_box(space_group_settings(19)));
    });
    group.finish();
}

fn bench_crystal_workflows(c: &mut Criterion) {
    let context = ExecutionContext::default();
    let Some(bytes) = Sample::Tiny.cif() else {
        panic!("tiny benchmark fixture has no mmCIF");
    };
    let buffer = input(bytes);
    let document = match molframe_cif::parse(&buffer) {
        Ok((document, _)) => document,
        Err(findings) => panic!("crystal benchmark parse failed: {findings:?}"),
    };
    let symmetry = match lower_symmetry(&document) {
        Ok(symmetry) => symmetry,
        Err(findings) => panic!("crystal benchmark symmetry failed: {findings:?}"),
    };
    let structure = structure(Sample::Tiny);
    if let Err(finding) = collect_crystal_neighbors(
        &structure,
        &symmetry,
        ModelIndex::new(0),
        5.0,
        molframe_xtal::CrystalNeighborOptions::default(),
        &context,
    ) {
        panic!("crystal-neighbour benchmark setup failed: {finding}");
    }
    let mut group = c.benchmark_group("xtal_workflows");
    group.bench_function("lower_symmetry", |b| {
        b.iter(|| black_box(lower_symmetry(&document)));
    });
    group.bench_function("lower_assemblies", |b| {
        b.iter(|| black_box(lower_assemblies(&document)));
    });
    group.bench_function("crystal_neighbors", |b| {
        b.iter(|| {
            black_box(collect_crystal_neighbors(
                &structure,
                &symmetry,
                ModelIndex::new(0),
                5.0,
                molframe_xtal::CrystalNeighborOptions::default(),
                &context,
            ))
        });
    });
    group.finish();
}

fn bench_extended_crystal_workflows(c: &mut Criterion) {
    let map = density_map();
    let map_bytes = map.to_mrc_bytes().required("MRC benchmark fixture failed");
    let mut block_reader = MrcBlockReader::new(
        Cursor::new(map_bytes.as_slice()),
        MrcBlockOptions::default(),
    )
    .required("MRC block-reader fixture failed");
    let mut block_values = Vec::new();
    let reflection = reflection_table();
    let mtz = write_mtz(&reflection).required("MTZ benchmark fixture failed");
    let structure_factor_document = parse_document(REFLECTION_CIF);
    let ncs_document = parse_document(NCS_CIF);
    let ncs = lower_ncs(&ncs_document).required("NCS benchmark fixture failed");
    let cell = UnitCell {
        lengths: [10.0, 11.0, 12.0],
        angles: [90.0, 91.0, 92.0],
    };
    let cell_transform = CellTransform::new(&cell).required("cell transform fixture failed");
    let affine = AffineTransform::new(
        [[1.0, 0.1, 0.0], [0.0, 1.0, 0.2], [0.0, 0.0, 1.0]],
        [1.0, 2.0, 3.0],
    );
    let mask: Vec<_> = (0..map.values.len()).map(|index| index % 3 != 0).collect();
    let mut group = c.benchmark_group("xtal_extended");
    group.bench_function("mrc_write/4096", |b| {
        b.iter(|| black_box(map.to_mrc_bytes()));
    });
    group.bench_function("mrc_read/4096", |b| {
        b.iter(|| black_box(DensityMap::from_mrc_bytes(&map_bytes)));
    });
    group.bench_function("mrc_block_read/8x8x8", |b| {
        b.iter(|| {
            block_reader
                .read_block_into([3, 3, 3], [8, 8, 8], &mut block_values)
                .required("MRC benchmark block failed");
            black_box(&block_values);
        });
    });
    group.bench_function("map_statistics/4096", |b| {
        b.iter(|| black_box(map.statistics()));
    });
    group.bench_function("map_masked_statistics/4096", |b| {
        b.iter(|| black_box(map.masked_statistics(&mask)));
    });
    group.bench_function("map_histogram/64x4096", |b| {
        b.iter(|| black_box(map.histogram(64, 0.0, 10.0)));
    });
    group.bench_function("cell_transform/1024", |b| {
        b.iter(|| {
            let fractional = [0.2, 0.3, 0.4];
            black_box(cell_transform.to_fractional(cell_transform.to_cartesian(fractional)))
        });
    });
    group.bench_function("affine_transform/1024", |b| {
        let points = [[0.0_f32, 1.0, 2.0]; 1_024];
        b.iter(|| {
            let next = affine.then(&AffineTransform::IDENTITY);
            black_box(
                points
                    .iter()
                    .map(|&point| next.apply(point))
                    .collect::<Vec<_>>(),
            )
        });
    });
    bench_reflection_workflows(
        &mut group,
        &ReflectionFixtures {
            reflection: &reflection,
            mtz: &mtz,
            structure_factor_document: &structure_factor_document,
            ncs_document: &ncs_document,
            ncs: &ncs,
            map: &map,
        },
    );
    group.finish();
}

fn bench_text_grids(c: &mut Criterion) {
    let [cube, dx] = text_grid_fixtures(64);
    let mut group = c.benchmark_group("xtal_text_grid_read");
    group.bench_function("cube/64x64x64", |b| {
        b.iter(|| black_box(read_cube(cube.as_bytes())));
    });
    group.bench_function("dx/64x64x64", |b| {
        b.iter(|| black_box(read_dx(dx.as_bytes())));
    });
    group.finish();
}

fn text_grid_fixtures(side: usize) -> [String; 2] {
    let values = side * side * side;
    let mut cube =
        format!("density\nbenchmark\n0 0 0 0\n{side} 1 0 0\n{side} 0 1 0\n{side} 0 0 1\n");
    let mut dx = format!(
        "object 1 class gridpositions counts {side} {side} {side}\n\
         origin 0 0 0\n\
         delta 1 0 0\n\
         delta 0 1 0\n\
         delta 0 0 1\n\
         object 2 class gridconnections counts {side} {side} {side}\n\
         object 3 class array type float rank 0 items {values} data follows\n"
    );
    for index in 0..values {
        let value = index % 97;
        write!(cube, "{value}.0 ").required("cube benchmark fixture formatting failed");
        write!(dx, "{value}.0 ").required("DX benchmark fixture formatting failed");
        if index % 8 == 7 {
            cube.push('\n');
            dx.push('\n');
        }
    }
    [cube, dx]
}

/// Everything the reflection and symmetry benchmarks read, built once by the
/// caller so no fixture is reconstructed inside a measured iteration.
struct ReflectionFixtures<'a> {
    reflection: &'a ReflectionTable,
    mtz: &'a [u8],
    structure_factor_document: &'a molframe_cif::Document,
    ncs_document: &'a molframe_cif::Document,
    ncs: &'a NcsSet,
    map: &'a DensityMap,
}

/// Reflection tables, MTZ round trips and NCS lowering, kept out of the density
/// benchmark so neither function outgrows a readable screenful.
fn bench_reflection_workflows(
    group: &mut BenchmarkGroup<'_, WallTime>,
    fixtures: &ReflectionFixtures<'_>,
) {
    let ReflectionFixtures {
        reflection,
        mtz,
        structure_factor_document,
        ncs_document,
        ncs,
        map,
    } = *fixtures;
    group.bench_function("reflection_validate/512", |b| {
        b.iter(|| black_box(reflection.validate()));
    });
    group.bench_function("reflection_miller_indices/512", |b| {
        b.iter(|| black_box(reflection.miller_indices()));
    });
    group.bench_function("mtz_write/512x4", |b| {
        b.iter(|| black_box(write_mtz(reflection)));
    });
    group.bench_function("mtz_read/512x4", |b| {
        b.iter(|| black_box(read_mtz(mtz)));
    });
    group.bench_function("structure_factor_lower/2", |b| {
        b.iter(|| black_box(lower_structure_factor_cif(structure_factor_document)));
    });
    group.bench_function("structure_factor_write/2x6", |b| {
        let table = lower_structure_factor_cif(structure_factor_document)
            .required("structure-factor fixture failed");
        b.iter(|| black_box(write_structure_factor_cif(&table)));
    });
    group.bench_function("ncs_lower/2", |b| {
        b.iter(|| black_box(lower_ncs(ncs_document)));
    });
    group.bench_function("ncs_iterators/2", |b| {
        b.iter(|| black_box((ncs.len(), ncs.operators().count(), ncs.generators().count())));
    });
    group.bench_function("map_sample/1024", |b| {
        let Some(sampler) = map.sampler() else {
            panic!("density map cell must be valid");
        };
        b.iter(|| {
            black_box(
                (0..1_024)
                    .map(|_| sampler.sample_cartesian([1.0, 2.0, 3.0], MapBoundary::Missing))
                    .collect::<Vec<_>>(),
            )
        });
    });
}

const NCS_CIF: &str = r"data_ncs
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
given-op given 'already deposited' 1 0 0 0 0 1 0 0 0 0 1 0
new-op generate 'missing copy' 1 0 0 5 0 1 0 6 0 0 1 7
";

const REFLECTION_CIF: &str = "data_sf\n_cell.length_a 10\n_cell.length_b 11\n_cell.length_c 12\n_cell.angle_alpha 90\n_cell.angle_beta 91\n_cell.angle_gamma 92\n_symmetry.Int_Tables_number 4\n_symmetry.space_group_name_H-M 'P 1 21 1'\nloop_\n_refln.index_h\n_refln.index_k\n_refln.index_l\n_refln.F_meas_au\n_refln.F_meas_sigma_au\n_refln.status\n0 0 1 10.5 0.2 o\n1 0 0 ? . f\n";

fn parse_document(text: &str) -> molframe_cif::Document {
    molframe_cif::parse(&input(text.as_bytes()))
        .map(|(document, _)| document)
        .required("xtal CIF fixture failed")
}

fn density_map() -> DensityMap {
    let values: Vec<f32> = (0_u16..4_096)
        .map(|index| f32::from(index % 97) * 0.1)
        .collect();
    DensityMap {
        dimensions: [16, 16, 16],
        starts: [0, 0, 0],
        sampling: [16, 16, 16],
        cell: UnitCell {
            lengths: [16.0, 16.0, 16.0],
            angles: [90.0, 90.0, 90.0],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values,
    }
}

fn reflection_table() -> ReflectionTable {
    let indices = |axis: u16| {
        (0_u16..512)
            .map(|row| ReflectionValue::Integer(i64::from((row / 17 + axis) % 13) - 6))
            .collect()
    };
    let reals = |scale: f64| {
        (0_u16..512)
            .map(|row| ReflectionValue::Real(scale * f64::from(row % 101)))
            .collect()
    };
    ReflectionTable {
        title: "bench".into(),
        cell: Some(UnitCell {
            lengths: [10.0, 11.0, 12.0],
            angles: [90.0, 91.0, 92.0],
        }),
        space_group_number: Some(4),
        space_group_name: Some("P 1 21 1".into()),
        columns: vec![
            ReflectionColumn {
                label: "H".into(),
                column_type: ReflectionColumnType::MillerIndex,
                values: indices(0),
                dataset_id: 0,
                mtz_type: Some('H'),
            },
            ReflectionColumn {
                label: "K".into(),
                column_type: ReflectionColumnType::MillerIndex,
                values: indices(1),
                dataset_id: 0,
                mtz_type: Some('H'),
            },
            ReflectionColumn {
                label: "L".into(),
                column_type: ReflectionColumnType::MillerIndex,
                values: indices(2),
                dataset_id: 0,
                mtz_type: Some('H'),
            },
            ReflectionColumn {
                label: "FP".into(),
                column_type: ReflectionColumnType::Amplitude,
                values: reals(0.5),
                dataset_id: 0,
                mtz_type: Some('F'),
            },
        ],
        datasets: Vec::new(),
        history: Vec::new(),
        symmetry_operations: vec!["X,Y,Z".into()],
        sort_order: [0; 5],
        resolution_range: None,
        missing_value: None,
        extra_header_records: Vec::new(),
    }
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_space_groups(&mut criterion);
    bench_crystal_workflows(&mut criterion);
    bench_extended_crystal_workflows(&mut criterion);
    bench_text_grids(&mut criterion);
    criterion.final_summary();
}
