//! Criterion coverage for the all-features facade read/write path.

use criterion::{Criterion, Throughput, black_box};
use molframe::{ReadOptions, read_bytes, write_mmcif};
use molframe_bench::{Sample, input};

fn bench_facade(c: &mut Criterion) {
    let bytes = match Sample::Medium.cif() {
        Some(bytes) => bytes.to_vec(),
        None => panic!("facade benchmark requires the medium CIF fixture"),
    };
    let options = ReadOptions::new();
    let input = input(&bytes);
    let (structure, _) = match molframe_cif::read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("facade fixture failed: {findings:?}"),
    };
    let mut group = c.benchmark_group("facade_io");
    group.throughput(Throughput::Bytes(bytes.len() as u64));
    group.bench_function("read_bytes", |b| {
        b.iter(|| black_box(read_bytes(bytes.clone(), Some("4hhb.cif"), &options)));
    });
    group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
    group.bench_function("write_mmcif", |b| {
        b.iter(|| black_box(write_mmcif(&structure)));
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_facade(&mut criterion);
    criterion.final_summary();
}
