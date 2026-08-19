//! Criterion coverage for `ModelCIF` lowering and canonical writing.

use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, input, structure};
use pdbiox_cif::parse;
use pdbiox_modelcif::{lower, write_canonical};

fn bench_lower(c: &mut Criterion) {
    let Some(bytes) = Sample::Tiny.cif() else {
        panic!("tiny fixture has no CIF bytes")
    };
    let buffer = input(bytes);
    let document = match parse(&buffer) {
        Ok((document, _)) => document,
        Err(findings) => panic!("ModelCIF bench fixture failed: {findings:?}"),
    };
    c.bench_function("modelcif_lower/1crn", |b| {
        b.iter(|| black_box(lower(&document)));
    });
}

fn bench_write(c: &mut Criterion) {
    let structure = structure(Sample::Tiny);
    let model_cif = pdbiox_modelcif::ModelCif::default();
    c.bench_function("modelcif_write/1crn", |b| {
        b.iter(|| black_box(write_canonical(&structure, &model_cif)));
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_lower(&mut criterion);
    bench_write(&mut criterion);
    criterion.final_summary();
}
