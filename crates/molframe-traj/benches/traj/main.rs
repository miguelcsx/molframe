//! Criterion coverage for trajectory parsing and ensemble operations.

use criterion::{BatchSize, Criterion, black_box};
use molframe_bench::{Sample, coordinates, perturbed, structure};
use molframe_traj::{
    CartesianFit, DcdReader, DcdWriteOptions, DielectricOptions, EnsembleDistanceMatrix,
    FrameAlignment, H5mdOptions, HarmonicSimilarityOptions, KMeansOptions, Linkage, NamdEndian,
    PathFrameMetric, RemainderPolicy, SurvivalMode, Timestep, TrrWriteOptions, XtcWriteOptions,
    agglomerative_clustering, block_convergence, cartesian_pca, cluster_population_similarity,
    dbscan_clustering, dielectric_from_dipoles, diffusion_map, generalized_procrustes_mean,
    group_coordinate_variance, harmonic_ensemble_similarity, kmeans, mean_squared_displacement,
    medoid, pairwise_fitted_rmsd, parse_aims_geometry, parse_charmm_record, parse_dcd,
    parse_gro_records, parse_gromacs_itp, parse_h5md, parse_hoomd_xml, parse_lammps_data,
    parse_lammps_dump, parse_namd_binary, parse_psf, parse_trr, parse_txyz_records, parse_xtc,
    parse_xyz, path_similarity, rmsd_to_reference, water_dynamics, write_aims_geometry,
    write_charmm_card, write_dcd, write_gro, write_h5md, write_namd_binary, write_psf, write_trr,
    write_txyz, write_xtc, write_xyz,
};

use std::{fmt::Debug, mem::size_of};

mod decode;

pub(crate) trait BenchRequired<T> {
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

const XYZ: &str =
    "3\nframe one\nC 0 0 0\nN 1 0 0\nO 0 1 0\n3\nframe two\nC 0.1 0 0\nN 1.1 0 0\nO 0 1.1 0\n";

const AIMS: &str = "lattice_vector 10 0 0\nlattice_vector 0 11 0\nlattice_vector 0 0 12\natom 0 0 0 C\natom 1 2 3 H\n";

const CHARMM: &str =
    "* title\n*\n    1\n    1    1 MOL  CA  -123.45600   2.50000  99.00000 SYS  1      0.00000\n";

const GRO: &str = "water\n2\n    1WAT     OW    1   1.000   2.000   3.000\n    1WAT     HW    2   1.100   2.000   3.000\n   5.00000   5.00000   5.00000\n";

const ITP: &str = "[ moleculetype ]\nSOL 2\n\n[ atoms ]\n1 OW 1 SOL OW 1 -0.834 15.9994\n2 HW 1 SOL HW 1 0.417 1.008\n\n[ bonds ]\n1 2 1 0.09572 345000\n";

const LAMMPS_DUMP: &str = "ITEM: TIMESTEP\n20\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 20\n0 30\nITEM: ATOMS id type x y z\n2 1 2 3 4\n1 1 1 2 3\n";

const LAMMPS_DATA: &str = "water topology\n\n2 atoms\n1 bonds\n\n0 10 xlo xhi\n0 11 ylo yhi\n0 12 zlo zhi\n\nMasses\n\n1 15.9994\n2 1.008\n\nAtoms # full\n\n1 1 1 -0.834 0 0 0 0 0 0\n2 1 2 0.417 1 0 0 0 0 0\n\nBonds\n\n1 1 1 2\n";

const PSF: &str = "PSF EXT XPLOR\n\n       1 !NTITLE\n REMARKS test\n\n       2 !NATOM\n       1 SEG 1 WAT O  O   -0.30 15.999 0\n       2 SEG 1 WAT H  H    0.30  1.008 0\n\n       1 !NBOND: bonds\n       1 2\n";

const TXYZ: &str = "2  water\n1 O 0.0 0.0 0.0 1 2\n2 H 1.0 0.0 0.0 2 1\n";

const HOOMD_XML: &str = "<?xml version=\"1.0\"?><hoomd_xml version=\"1.6\"><configuration time_step=\"42\" dimensions=\"3\" natoms=\"2\"><box lx=\"10\" ly=\"11\" lz=\"12\"/><position>0 0 0 1 0 0</position><type>O H</type></configuration></hoomd_xml>";

fn bench_formats(c: &mut Criterion) {
    let Some(frames) = parse_xyz(XYZ) else {
        panic!("embedded XYZ benchmark fixture is invalid")
    };
    let mut group = c.benchmark_group("traj_formats");
    group.bench_function("parse_xyz", |b| b.iter(|| black_box(parse_xyz(XYZ))));
    group.bench_function("write_xyz", |b| b.iter(|| black_box(write_xyz(&frames))));
    group.finish();
}

fn bench_ensemble(c: &mut Criterion) {
    let structure = structure(Sample::Tiny);
    let reference = coordinates(&structure);
    let model = perturbed(&reference, 0.01);
    let frames = vec![
        Timestep {
            frame: 0,
            positions: reference.clone(),
            ..Timestep::default()
        },
        Timestep {
            frame: 1,
            positions: model,
            ..Timestep::default()
        },
    ];
    c.bench_function("traj_rmsd_to_reference", |b| {
        b.iter(|| black_box(rmsd_to_reference(&frames, 0, FrameAlignment::Rigid)));
    });

    let coordinate_frames = synthetic_frames(&reference[..reference.len().min(32)], 16);
    let streaming_path_frames = synthetic_frames(&[[0.0, 0.0, 0.0]], 1_024);
    let streaming_path_workspace = 2 * 512 * size_of::<f64>();
    let mut group = c.benchmark_group("traj_ensemble_kernels");
    group.bench_function("mean_squared_displacement", |b| {
        b.iter(|| black_box(mean_squared_displacement(&coordinate_frames, &[], 8)));
    });
    group.bench_function("cartesian_pca", |b| {
        b.iter(|| {
            black_box(cartesian_pca(
                &coordinate_frames,
                CartesianFit::None,
                3,
                32 * 1_024 * 1_024,
            ))
        });
    });
    group.bench_function("path_similarity/cartesian", |b| {
        b.iter(|| {
            black_box(path_similarity(
                &coordinate_frames[..8],
                &coordinate_frames[8..],
                PathFrameMetric::CartesianRmsd,
                1_024 * 1_024,
            ))
        });
    });
    group.bench_function("path_similarity/fitted", |b| {
        b.iter(|| {
            black_box(path_similarity(
                &coordinate_frames[..8],
                &coordinate_frames[8..],
                PathFrameMetric::FittedRmsd,
                1_024 * 1_024,
            ))
        });
    });
    group.bench_function("path_similarity/cartesian_streaming_512x512x1", |b| {
        b.iter(|| {
            black_box(
                path_similarity(
                    &streaming_path_frames[..512],
                    &streaming_path_frames[512..],
                    PathFrameMetric::CartesianRmsd,
                    streaming_path_workspace,
                )
                .required("streaming path workspace exceeded its linear ceiling"),
            )
        });
    });
    group.finish();
}

fn bench_clustering(c: &mut Criterion) {
    let observations: Vec<Vec<f64>> = (0..2_048)
        .map(|row| (0..8).map(|column| synthetic_value(row, column)).collect())
        .collect();
    let options = KMeansOptions {
        initial_centres: &[0, 1_024],
        maximum_iterations: 50,
        convergence_tolerance_squared: 1e-12,
    };
    c.bench_function("traj_clustering/kmeans_2048x8", |b| {
        b.iter(|| black_box(kmeans(&observations, options)));
    });

    let size = 128_usize;
    let values: Vec<f64> = (0..size)
        .flat_map(|left| {
            (0..size).map(move |right| {
                let delta = left.abs_diff(right);
                f64::from(usize_to_u32(delta))
            })
        })
        .collect();
    let distances = EnsembleDistanceMatrix {
        size,
        values: values.into_boxed_slice(),
    };
    let members: Vec<_> = (0..size).collect();
    let mut group = c.benchmark_group("traj_clustering");
    group.bench_function("agglomerative_128", |b| {
        b.iter(|| black_box(agglomerative_clustering(&distances, 8, Linkage::Average)));
    });
    group.bench_function("dbscan_128", |b| {
        b.iter(|| black_box(dbscan_clustering(&distances, 4.0, 3)));
    });
    group.bench_function("medoid_128", |b| {
        b.iter(|| black_box(medoid(&distances, &members)));
    });
    group.finish();
}

fn bench_extended_ensemble(c: &mut Criterion) {
    let structure = structure(Sample::Tiny);
    let reference = coordinates(&structure);
    let frames = synthetic_frames(&reference[..reference.len().min(8)], 16);
    let distances = distance_matrix(32);
    let observations: Vec<Vec<f64>> = (0_usize..64)
        .map(|row| {
            (0..4)
                .map(|column| f64::from(usize_to_u32((row + 3 * column) % 17)) / 17.0)
                .collect()
        })
        .collect();
    let dipoles: Vec<[f64; 3]> = (0_usize..1_024)
        .map(|index| {
            let value = f64::from(usize_to_u32(index % 23)) / 23.0;
            [value, value * 0.5, -value * 0.25]
        })
        .collect();
    let groups = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];
    let values: Vec<f64> = (0_u32..1_024).map(f64::from).collect();
    let occupancy: Vec<Vec<bool>> = (0..1_024)
        .map(|frame| (0..16).map(|site| (frame + site) % 7 != 0).collect())
        .collect();
    let mut group = c.benchmark_group("traj_extended_ensemble");
    group.bench_function("pairwise_fitted_rmsd/16x8", |b| {
        b.iter(|| black_box(pairwise_fitted_rmsd(&frames, 32 * 1_024 * 1_024)));
    });
    group.bench_function("generalized_procrustes/16x8", |b| {
        b.iter(|| black_box(generalized_procrustes_mean(&frames, 1.0e-6, 10)));
    });
    group.bench_function("diffusion_map/32x4", |b| {
        b.iter(|| black_box(diffusion_map(&distances, 2.0, 2, 4)));
    });
    group.bench_function("harmonic_similarity/64x4", |b| {
        b.iter(|| {
            black_box(harmonic_ensemble_similarity(
                &observations,
                &observations,
                HarmonicSimilarityOptions {
                    covariance_regularization: 1.0e-6,
                    memory_limit_bytes: 1_024 * 1_024,
                },
            ))
        });
    });
    group.bench_function("group_coordinate_variance/16x8", |b| {
        b.iter(|| black_box(group_coordinate_variance(&frames, &groups)));
    });
    group.bench_function("block_convergence/1024", |b| {
        b.iter(|| black_box(block_convergence(&values, 32, RemainderPolicy::Reject)));
    });
    group.bench_function("dielectric/1024", |b| {
        b.iter(|| {
            black_box(dielectric_from_dipoles(
                &dipoles,
                DielectricOptions {
                    fluctuation_prefactor: 2.0,
                },
            ))
        });
    });
    group.bench_function("water_dynamics/1024x16", |b| {
        b.iter(|| {
            black_box(water_dynamics(
                &occupancy,
                16,
                0.5,
                SurvivalMode::Continuous,
            ))
        });
    });
    group.bench_function("cluster_population_similarity/1024", |b| {
        let labels: Vec<_> = (0..1_024).map(|value| value % 8).collect();
        b.iter(|| black_box(cluster_population_similarity(&labels, &labels, 8)));
    });
    group.finish();
}

fn bench_extended_formats(c: &mut Criterion) {
    let aims = parse_aims_geometry(AIMS).required("AIMS fixture failed");
    let charmm = parse_charmm_record(CHARMM).required("CHARMM fixture failed");
    let gro = parse_gro_records(GRO).required("GRO fixture failed");
    let _itp = parse_gromacs_itp(ITP).required("ITP fixture failed");
    let _lammps_dump = parse_lammps_dump(LAMMPS_DUMP).required("LAMMPS dump fixture failed");
    let _lammps_data = parse_lammps_data(LAMMPS_DATA).required("LAMMPS data fixture failed");
    let psf = parse_psf(PSF).required("PSF fixture failed");
    let txyz = parse_txyz_records(TXYZ).required("TXYZ fixture failed");
    let _hoomd = parse_hoomd_xml(HOOMD_XML).required("HOOMD fixture failed");
    let binary_frames: Vec<_> = synthetic_frames(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        4,
    )
    .into_iter()
    .enumerate()
    .map(|(frame, positions)| Timestep {
        frame,
        time: Some(f64::from(usize_to_u32(frame))),
        positions,
        ..Timestep::default()
    })
    .collect();
    let dcd = write_dcd(&binary_frames, &DcdWriteOptions::default()).required("DCD fixture failed");
    let trr = write_trr(&binary_frames, TrrWriteOptions::default()).required("TRR fixture failed");
    let xtc = write_xtc(&binary_frames, XtcWriteOptions::default()).required("XTC fixture failed");
    let namd =
        write_namd_binary(&binary_frames[0], NamdEndian::Little).required("NAMD fixture failed");
    let h5md = write_h5md(&binary_frames, &H5mdOptions::default()).required("H5MD fixture failed");
    let mut group = c.benchmark_group("traj_extended_formats");
    group.bench_function("parse/aims", |b| {
        b.iter(|| black_box(parse_aims_geometry(AIMS)));
    });
    group.bench_function("write/aims", |b| {
        b.iter(|| black_box(write_aims_geometry(&aims)));
    });
    group.bench_function("parse/charmm", |b| {
        b.iter(|| black_box(parse_charmm_record(CHARMM)));
    });
    group.bench_function("write/charmm", |b| {
        b.iter(|| black_box(write_charmm_card(&charmm)));
    });
    group.bench_function("parse/gro", |b| {
        b.iter(|| black_box(parse_gro_records(GRO)));
    });
    group.bench_function("write/gro", |b| b.iter(|| black_box(write_gro(&gro))));
    group.bench_function("parse/gromacs_itp", |b| {
        b.iter(|| black_box(parse_gromacs_itp(ITP)));
    });
    group.bench_function("parse/lammps_dump", |b| {
        b.iter(|| black_box(parse_lammps_dump(LAMMPS_DUMP)));
    });
    group.bench_function("parse/lammps_data", |b| {
        b.iter(|| black_box(parse_lammps_data(LAMMPS_DATA)));
    });
    group.bench_function("parse/psf", |b| b.iter(|| black_box(parse_psf(PSF))));
    group.bench_function("write/psf", |b| b.iter(|| black_box(write_psf(&psf))));
    group.bench_function("parse/txyz", |b| {
        b.iter(|| black_box(parse_txyz_records(TXYZ)));
    });
    group.bench_function("write/txyz", |b| b.iter(|| black_box(write_txyz(&txyz))));
    group.bench_function("parse/hoomd_xml", |b| {
        b.iter(|| black_box(parse_hoomd_xml(HOOMD_XML)));
    });
    group.bench_function("write/dcd/4x4", |b| {
        b.iter(|| black_box(write_dcd(&binary_frames, &DcdWriteOptions::default())));
    });
    group.bench_function("parse/dcd/4x4", |b| b.iter(|| black_box(parse_dcd(&dcd))));
    group.bench_function("write/trr/4x4", |b| {
        b.iter(|| black_box(write_trr(&binary_frames, TrrWriteOptions::default())));
    });
    group.bench_function("parse/trr/4x4", |b| b.iter(|| black_box(parse_trr(&trr))));
    group.bench_function("write/xtc/4x4", |b| {
        b.iter(|| black_box(write_xtc(&binary_frames, XtcWriteOptions::default())));
    });
    group.bench_function("parse/xtc/4x4", |b| b.iter(|| black_box(parse_xtc(&xtc))));
    group.bench_function("write/namd/4", |b| {
        b.iter(|| black_box(write_namd_binary(&binary_frames[0], NamdEndian::Little)));
    });
    group.bench_function("parse/namd/4", |b| {
        b.iter(|| black_box(parse_namd_binary(&namd)));
    });
    group.bench_function("write/h5md/4x4", |b| {
        b.iter(|| black_box(write_h5md(&binary_frames, &H5mdOptions::default())));
    });
    group.bench_function("parse/h5md/4x4", |b| {
        b.iter(|| black_box(parse_h5md(&h5md)));
    });
    group.finish();
}

fn bench_dcd_streaming(c: &mut Criterion) {
    let sample = structure(Sample::Small);
    let frames = synthetic_frames(&coordinates(&sample), 256)
        .into_iter()
        .enumerate()
        .map(|(frame, positions)| Timestep {
            frame,
            positions,
            ..Timestep::default()
        })
        .collect::<Vec<_>>();
    let bytes = write_dcd(&frames, &DcdWriteOptions::default()).required("DCD fixture failed");
    let file = tempfile::NamedTempFile::new().required("DCD temporary file failed");
    std::fs::write(file.path(), &bytes).required("DCD fixture write failed");
    let path = file.path().to_path_buf();
    let mut group = c.benchmark_group("traj_dcd_bounded");
    group.sample_size(60);
    group.bench_function("materialized/256_frames", |b| {
        b.iter(|| {
            let input = std::fs::read(&path).required("DCD benchmark read failed");
            black_box(parse_dcd(&input).required("DCD materialized parse failed"));
        });
    });
    group.bench_function("streaming/256_frames", |b| {
        b.iter_batched(
            || DcdReader::open(&path).required("DCD streaming open failed"),
            |mut reader| {
                let mut frame = Timestep::default();
                let mut count = 0_usize;
                while reader
                    .read_next_frame(&mut frame)
                    .required("DCD streaming read failed")
                {
                    count += 1;
                    black_box(&frame.positions);
                }
                black_box(count);
            },
            BatchSize::SmallInput,
        );
    });
    group.finish();
}

fn distance_matrix(size: usize) -> EnsembleDistanceMatrix {
    let values: Vec<f64> = (0..size)
        .flat_map(|left| {
            (0..size).map(move |right| {
                let delta = left.abs_diff(right);
                f64::from(usize_to_u32(delta))
            })
        })
        .collect();
    EnsembleDistanceMatrix {
        size,
        values: values.into_boxed_slice(),
    }
}

pub(crate) fn synthetic_frames(reference: &[[f32; 3]], frame_count: usize) -> Vec<Vec<[f32; 3]>> {
    (0..frame_count)
        .map(|frame| {
            let phase = f32::from(usize_to_u16(frame)) * 0.002;
            reference
                .iter()
                .enumerate()
                .map(|(atom, &[x, y, z])| {
                    let atom_phase = f32::from(usize_to_u16(atom % 17));
                    [x + phase, y - phase * 0.5, z + phase * atom_phase * 0.01]
                })
                .collect()
        })
        .collect()
}

fn synthetic_value(row: usize, column: usize) -> f64 {
    let cluster = if row < 1_024 { 0.0 } else { 10.0 };
    let residue = (row * 31 + column * 17) % 101;
    cluster + f64::from(usize_to_u32(residue)) / 101.0
}

pub(crate) fn usize_to_u32(value: usize) -> u32 {
    let Ok(value) = u32::try_from(value) else {
        return u32::MAX;
    };
    value
}

fn usize_to_u16(value: usize) -> u16 {
    let Ok(value) = u16::try_from(value) else {
        return u16::MAX;
    };
    value
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_formats(&mut criterion);
    bench_ensemble(&mut criterion);
    bench_clustering(&mut criterion);
    bench_extended_ensemble(&mut criterion);
    bench_extended_formats(&mut criterion);
    bench_dcd_streaming(&mut criterion);
    decode::bench_decode_throughput(&mut criterion);
    criterion.final_summary();
}
