//! Benchmarks for `molframe-pdb` hot paths.
//!
//! Covers legacy PDB reading and writing, MMTF round-tripping, and the hybrid-36
//! integer codec used for large serials and residue numbers.

use std::fmt::Write;

use criterion::{Criterion, Throughput, black_box};
use molframe_bench::{Sample, input, structure, structure_from_pdb};
use molframe_core::io::ReadOptions;
use molframe_pdb::{PdbOptions, hybrid36, read, write, write_mmtf};

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

/// Fixed-column PDB text of `atoms` atoms: four atoms per residue, 9999
/// residues per chain, so no serial, residue or chain label overflows.
fn synthetic_pdb(atoms: u32) -> String {
    let mut text = String::new();
    for index in 0..atoms {
        let name = ["N", "CA", "C", "O"][(index % 4) as usize];
        let residue = index / 4;
        let Some(chain) = char::from_u32(u32::from(b'A') + (residue / 9999) % 26) else {
            break;
        };
        let x = f64::from(index % 997) + 0.5;
        let y = f64::from(index % 991) * 0.25;
        let z = f64::from(index % 983) * -0.125;
        let serial = index % 99_999 + 1;
        let seq = residue % 9999 + 1;
        let _ = writeln!(
            text,
            "ATOM  {serial:>5}  {name:<3} GLY {chain}{seq:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00           {}",
            &name[..1],
        );
    }
    text.push_str("END\n");
    text
}

/// Read time against atom count, to keep the per-atom cost flat at 1M atoms.
fn bench_read_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("pdb_read_atoms");
    group.sample_size(10);
    for atoms in [100_000_u32, 1_000_000] {
        let text = synthetic_pdb(atoms);
        let input = input(text.as_bytes());
        group.throughput(Throughput::Elements(u64::from(atoms)));
        group.bench_function(format!("{atoms}"), |b| {
            b.iter(|| {
                let _ = black_box(read(&input, &ReadOptions::new()));
            });
        });
    }
    group.finish();
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
    bench_read_scaling(&mut criterion);
    bench_write(&mut criterion);
    bench_mmtf_refusal(&mut criterion);
    bench_hybrid36(&mut criterion);
    criterion.final_summary();
}
