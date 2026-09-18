//! Benchmarks for `molframe-bcif` hot paths.
//!
//! Covers `BinaryCIF` reading, container-level decoding, serialization, and the
//! integer/float/string column encoders.

use criterion::{Criterion, Throughput, black_box};
use molframe_bcif::{
    BinaryDocument, DataType, EncodedData, Encoding, decode, encode_floats, encode_integers,
    encode_interval, encode_strings, read, write_structure,
};
use molframe_bench::{Sample, input, structure};
use molframe_core::io::{Limits, ReadOptions};

fn bench_read(c: &mut Criterion) {
    let mut group = c.benchmark_group("bcif_read");
    for sample in Sample::MODELS {
        let buffer = input(sample.bcif());
        group.throughput(Throughput::Bytes(sample.bcif().len() as u64));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(read(&buffer, &ReadOptions::new()));
            });
        });
    }
    group.finish();
}

fn bench_container_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("bcif_parse_container");
    for sample in Sample::MODELS {
        let bytes = sample.bcif();
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(BinaryDocument::parse(bytes, Limits::default()));
            });
        });
    }
    group.finish();
}

fn bench_column_decode(c: &mut Criterion) {
    // Build a realistic integer column by re-encoding the medium structure's
    // atom-site ids, then measure only the decode step.
    let structure = structure(Sample::Medium);
    let ids: Vec<i64> = (0..structure.atom_count()).map(i64::from).collect();
    let encoded = match encode_integers(&ids) {
        Ok(encoded) => encoded,
        Err(finding) => panic!("failed to encode integer column: {finding:?}"),
    };

    let mut group = c.benchmark_group("bcif_decode_integers");
    group.throughput(Throughput::Elements(ids.len() as u64));
    group.bench_function("4hhb", |b| {
        b.iter(|| {
            let _ = black_box(decode(&encoded));
        });
    });

    let wide_data: Vec<u8> = (0_i32..8_000_000).flat_map(i32::to_le_bytes).collect();
    let wide_encoded = EncodedData {
        encoding: vec![Encoding::ByteArray {
            r#type: DataType::Int32,
        }],
        data: wide_data,
    };
    group.throughput(Throughput::Elements(8_000_000));
    group.bench_function("wide_payload", |b| {
        b.iter(|| {
            let _ = black_box(decode(&wide_encoded));
        });
    });
    group.finish();
}

fn bench_write(c: &mut Criterion) {
    let mut group = c.benchmark_group("bcif_write_structure");
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium, Sample::Large] {
        let structure = structure(sample);
        group.throughput(Throughput::Elements(u64::from(structure.atom_count())));
        group.bench_function(sample.label(), |b| {
            b.iter(|| {
                let _ = black_box(write_structure(&structure));
            });
        });
    }
    group.finish();
}

fn bench_encoders(c: &mut Criterion) {
    let ints: Vec<i64> = (0..10_000).map(i64::from).collect();
    let floats: Vec<f64> = (0..10_000).map(|i| f64::from(i) * 0.1).collect();
    let strings: Vec<String> = (0..1_000).map(|i| format!("atom{i:04}")).collect();

    let mut group = c.benchmark_group("bcif_encode");
    group.throughput(Throughput::Elements(10_000));
    group.bench_function("integers", |b| {
        b.iter(|| {
            let _ = black_box(encode_integers(&ints));
        });
    });
    group.bench_function("floats", |b| {
        b.iter(|| {
            let _ = black_box(encode_floats(&floats));
        });
    });
    group.bench_function("strings", |b| {
        b.iter(|| {
            let _ = black_box(encode_strings(&strings));
        });
    });
    group.bench_function("interval", |b| {
        b.iter(|| {
            let _ = black_box(encode_interval(&floats, 100));
        });
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_read(&mut criterion);
    bench_container_parse(&mut criterion);
    bench_column_decode(&mut criterion);
    bench_write(&mut criterion);
    bench_encoders(&mut criterion);
    criterion.final_summary();
}
