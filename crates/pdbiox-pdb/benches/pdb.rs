//! Benchmarks for `pdbiox-pdb` hot paths.
//!
//! Covers legacy PDB reading and writing, MMTF round-tripping, and the hybrid-36
//! integer codec used for large serials and residue numbers.

use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, input, structure, structure_from_pdb};
use pdbiox_core::io::ReadOptions;
use pdbiox_pdb::{PdbOptions, hybrid36, read, write, write_mmtf};

fn bench_read(c: &mut Criterion) {
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let Some(bytes) = sample.pdb() else {
            panic!("sample {} has no PDB fixture", sample.label())
        };
        let input = input(bytes);
        c.bench_function(&format!("pdb_read/{}", sample.label()), |b| {
            b.iter(|| {
                let _ = black_box(read(&input, &ReadOptions::new()));
            });
        });
    }
}

fn bench_write(c: &mut Criterion) {
    let options = PdbOptions::new();
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let Some(structure) = structure_from_pdb(sample) else {
            panic!("sample {} has no PDB structure", sample.label())
        };
        c.bench_function(&format!("pdb_write/{}", sample.label()), |b| {
            b.iter(|| {
                let _ = black_box(write(&structure, &options));
            });
        });
    }
}

fn bench_mmtf_refusal(c: &mut Criterion) {
    let structure = structure(Sample::Tiny);
    c.bench_function("pdb_write_mmtf/refusal_without_metadata", |b| {
        b.iter(|| {
            let _ = black_box(write_mmtf(&structure));
        });
    });
}

fn bench_hybrid36(c: &mut Criterion) {
    let mut group = c.benchmark_group("pdb_hybrid36");
    group.bench_function("decode", |b| {
        b.iter(|| {
            for value in ["    A1", "999999", "   1A0", "   9ZZ"] {
                let _ = black_box(hybrid36::decode(value, 6));
            }
        });
    });
    group.bench_function("encode", |b| {
        b.iter(|| {
            for value in [1_i64, 99_999, 1_000_000, 2_176_782] {
                let _ = black_box(hybrid36::encode(value, 6));
            }
        });
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_read(&mut criterion);
    bench_write(&mut criterion);
    bench_mmtf_refusal(&mut criterion);
    bench_hybrid36(&mut criterion);
    criterion.final_summary();
}
