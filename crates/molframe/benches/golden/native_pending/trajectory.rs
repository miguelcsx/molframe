use std::fmt::Debug;

use criterion::{BenchmarkGroup, Throughput, black_box};
use molframe::{AnalysisPolicy, AtomSelection, SpatialBackend};
use molframe_bench::{Sample, structure};

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

const XYZ: &str = r"2
frame-0
C 0 0 0
N 1 0 0
2
frame-1
C 1 0 0
N 2 0 0
2
frame-2
C 2 0 0
N 3 0 0
";

pub(super) fn register(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    bench_gw_024(group);
    bench_gw_025(group);
    bench_gw_026(group);
}

fn bench_gw_024(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let positions = [
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.1, 0.0, 0.0], [1.1, 0.0, 0.0]],
    ];
    let selection = AtomSelection::All(2);
    let options = molframe::analysis::RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 2.0,
        bins: 4,
        volume: 100.0,
        backend: SpatialBackend::BruteForce,
    };
    group.throughput(Throughput::Elements(4));
    group.bench_function("GW-024", |b| {
        b.iter(|| {
            let mut neighbours =
                molframe::traj::FrameNeighborList::new([0_u32, 1], [0_u32, 1], 2.0, 0.5);
            let counts: Vec<_> = positions
                .into_iter()
                .map(|frame| {
                    neighbours
                        .pairs_positions(&frame, None)
                        .required("GW-024 neighbour query failed");
                    molframe::analysis::radial_distribution(
                        &frame,
                        &selection,
                        &selection,
                        options,
                        None,
                        &molframe::ExecutionContext::default(),
                    )
                    .required("GW-024 RDF failed")
                    .iter()
                    .map(|bin| bin.count)
                    .sum::<u64>()
                })
                .collect();
            black_box((counts, neighbours.statistics()));
        });
    });
}

fn bench_gw_025(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let frames = molframe::traj::parse_xyz(XYZ).required("GW-025 setup failed");
    let selected: Vec<_> = [0_usize, 2]
        .into_iter()
        .map(|index| frames[index].clone())
        .collect();
    group.throughput(Throughput::Bytes(XYZ.len() as u64));
    group.bench_function("GW-025", |b| {
        b.iter(|| {
            let rendered = molframe::traj::write_xyz(&selected);
            let round_trip =
                molframe::traj::parse_xyz(&rendered).required("GW-025 round trip failed");
            black_box(round_trip.len());
        });
    });
}

fn bench_gw_026(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = structure(Sample::Tiny);
    let trajectory = molframe::traj::Trajectory::from_frames(
        (0_u16..128)
            .map(|frame| molframe::traj::Frame {
                positions: structure
                    .positions()
                    .iter()
                    .map(|&[x, y, z]| [x + f32::from(frame) * 1.0e-4, y, z])
                    .collect(),
            })
            .collect(),
    )
    .required("GW-026 trajectory construction failed");
    let policy = AnalysisPolicy::default();
    let kernel = molframe::analysis::contacts_kernel(3.0, SpatialBackend::Auto);
    group.throughput(Throughput::Elements(128));
    group.bench_function("GW-026", |b| {
        b.iter(|| {
            let outputs: Vec<_> = [1_usize, 2, 4, 16]
                .into_iter()
                .map(|workers| {
                    let context = molframe::ExecutionContext::builder()
                        .worker_budget(workers)
                        .build()
                        .required("GW-026 execution context failed");
                    molframe::analysis::analyse_trajectory(
                        &structure,
                        &trajectory,
                        &policy,
                        &kernel,
                        &context,
                    )
                    .required("GW-026 failed")
                    .value
                })
                .collect();
            black_box(outputs);
        });
    });
}
