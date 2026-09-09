//! Criterion coverage for sequence alignment and sequence formats.

use criterion::{BenchmarkId, Criterion, Throughput, black_box};
use pdbiox_seq::{
    MsaOptions, Scoring, global, global_banded, kmer_counts, local, minimizers, parse_fasta,
    parse_fastq, progressive_msa, semi_global, write_fasta,
};

const FASTA: &str = ">a\nMKTAYIAKQRQISFVKSHFSRQ\n>b\nMKTAYIAKQ--ISFVKSHFSRQ\n";
const FASTQ: &str = "@a\nMKTAYIAKQRQ\n+\nIIIIIIIIIII\n";

fn bench_alignment(c: &mut Criterion) {
    let scoring = Scoring::simple();
    let mut group = c.benchmark_group("seq_alignment");
    for length in [64_usize, 256, 1_024] {
        let left = synthetic_sequence(length);
        let mut right = left.clone();
        for index in (17..right.len()).step_by(31) {
            right[index] = b'V';
        }
        group.throughput(Throughput::Elements(length as u64));
        group.bench_with_input(BenchmarkId::new("global", length), &length, |b, _| {
            b.iter(|| black_box(global(&left, &right, scoring)));
        });
        group.bench_with_input(BenchmarkId::new("local", length), &length, |b, _| {
            b.iter(|| black_box(local(&left, &right, scoring)));
        });
        group.bench_with_input(BenchmarkId::new("semi_global", length), &length, |b, _| {
            b.iter(|| black_box(semi_global(&left, &right, scoring)));
        });
        group.bench_with_input(BenchmarkId::new("banded", length), &length, |b, _| {
            b.iter(|| black_box(global_banded(&left, &right, scoring, 16)));
        });
    }
    group.finish();
}

fn bench_sequence_algorithms(c: &mut Criterion) {
    let sequence = synthetic_sequence(8_192);
    let mut group = c.benchmark_group("seq_indexing");
    group.throughput(Throughput::Bytes(sequence.len() as u64));
    group.bench_function("kmer_counts/k7", |b| {
        b.iter(|| black_box(kmer_counts(&sequence, 7)));
    });
    group.bench_function("minimizers/k7_w12", |b| {
        b.iter(|| black_box(minimizers(&sequence, 7, 12)));
    });
    group.finish();

    let sequences = [
        synthetic_sequence(96),
        synthetic_sequence(93),
        synthetic_sequence(89),
        synthetic_sequence(101),
    ];
    let sequence_views = sequences.each_ref().map(Vec::as_slice);
    c.bench_function("seq_msa/progressive_four", |b| {
        b.iter(|| {
            black_box(progressive_msa(
                &sequence_views,
                MsaOptions::progressive(Scoring::simple()),
            ))
        });
    });

    let larger: [Vec<u8>; 8] = std::array::from_fn(|index| {
        let mut sequence = synthetic_sequence(512 + index * 3);
        for position in (index + 11..sequence.len()).step_by(47 + index) {
            sequence[position] = b'Y';
        }
        sequence
    });
    let larger_views = larger.each_ref().map(Vec::as_slice);
    c.bench_function("seq_msa/progressive_eight_512", |b| {
        b.iter(|| {
            black_box(progressive_msa(
                &larger_views,
                MsaOptions::progressive(Scoring::simple()),
            ))
        });
    });

    let many = (0..64)
        .map(|index| {
            let mut sequence = synthetic_sequence(128);
            for position in (index % 23..sequence.len()).step_by(29 + index % 7) {
                sequence[position] = b'W';
            }
            sequence
        })
        .collect::<Vec<_>>();
    let many_views = many.iter().map(Vec::as_slice).collect::<Vec<_>>();
    c.bench_function("seq_msa/progressive_sixty_four_128", |b| {
        b.iter(|| {
            black_box(progressive_msa(
                &many_views,
                MsaOptions::progressive(Scoring::simple()),
            ))
        });
    });
}

fn bench_formats(c: &mut Criterion) {
    let mut group = c.benchmark_group("seq_formats");
    group.bench_function("fasta_round_trip", |b| {
        b.iter(|| {
            let records = parse_fasta(FASTA);
            black_box(write_fasta(&records))
        });
    });
    group.bench_function("fastq_parse", |b| b.iter(|| black_box(parse_fastq(FASTQ))));
    group.finish();
}

fn synthetic_sequence(length: usize) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    (0..length)
        .map(|index| ALPHABET[(index * 17 + index / 7) % ALPHABET.len()])
        .collect()
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_alignment(&mut criterion);
    bench_sequence_algorithms(&mut criterion);
    bench_formats(&mut criterion);
    criterion.final_summary();
}
