//! Benchmarks for `pdbiox-cif` hot paths.
//!
//! Covers the mmCIF read pipeline (parse, lower, full read, read with document
//! retention), canonical serialization, and PDBML decode.

use criterion::{Criterion, Throughput, black_box};
use pdbiox_bench::{Sample, input, structure_from_cif};
use pdbiox_cif::{lower, parse, read, read_pdbml, read_with_document, write_canonical};
use pdbiox_core::io::ReadOptions;

fn bench_read(c: &mut Criterion) {
    let mut group = c.benchmark_group("cif_read");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let Some(bytes) = sample.cif() else {
            panic!("sample {} has no CIF fixture", sample.label())
        };
        let input = input(bytes);
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(read(&input, &ReadOptions::new()));
            });
        });
    }
    group.finish();
}

fn bench_read_with_document(c: &mut Criterion) {
    let mut group = c.benchmark_group("cif_read_with_document");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let Some(bytes) = sample.cif() else {
            panic!("sample {} has no CIF fixture", sample.label())
        };
        let input = input(bytes);
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(read_with_document(&input, &ReadOptions::new()));
            });
        });
    }
    group.finish();
}

fn bench_tokenize_and_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("cif_parse");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let Some(bytes) = sample.cif() else {
            panic!("sample {} has no CIF fixture", sample.label())
        };
        let input = input(bytes);
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(parse(&input));
            });
        });
    }
    group.finish();
}

fn bench_lower(c: &mut Criterion) {
    let mut group = c.benchmark_group("cif_lower");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let Some(bytes) = sample.cif() else {
            panic!("sample {} has no CIF fixture", sample.label())
        };
        let input = input(bytes);
        let document = match parse(&input) {
            Ok((doc, _)) => doc,
            Err(findings) => panic!(
                "bench fixture {} parse failed: {findings:?}",
                sample.label()
            ),
        };
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(lower(&document, &ReadOptions::new()));
            });
        });
    }
    group.finish();
}

fn bench_write_canonical(c: &mut Criterion) {
    let mut group = c.benchmark_group("cif_write_canonical");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let Some(structure) = structure_from_cif(sample) else {
            panic!("sample {} has no CIF structure", sample.label())
        };
        let Some(bytes) = sample.cif() else {
            panic!("sample {} has no CIF fixture", sample.label())
        };
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(write_canonical(&structure));
            });
        });
    }
    group.finish();
}

fn bench_pdbml_read(c: &mut Criterion) {
    // PDBML is not committed; round-trip a tiny mmCIF through PDBML and measure
    // only the reader half.
    let Some(bytes) = Sample::Tiny.cif() else {
        panic!("tiny sample has no CIF fixture")
    };
    let input = input(bytes);
    let document = match parse(&input) {
        Ok((doc, _)) => doc,
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    };
    let xml = match pdbiox_cif::write_pdbml(&document) {
        Ok(bytes) => bytes,
        Err(finding) => panic!("failed to write PDBML for bench fixture: {finding:?}"),
    };
    let options = ReadOptions::new();
    let mut group = c.benchmark_group("cif_read_pdbml");
    group.throughput(Throughput::Bytes(xml.len() as u64));
    group.bench_function("1crn", |b| {
        b.iter(|| {
            let _ = black_box(read_pdbml(xml.as_bytes(), &options));
        });
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_read(&mut criterion);
    bench_read_with_document(&mut criterion);
    bench_tokenize_and_parse(&mut criterion);
    bench_lower(&mut criterion);
    bench_write_canonical(&mut criterion);
    bench_pdbml_read(&mut criterion);
    criterion.final_summary();
}
