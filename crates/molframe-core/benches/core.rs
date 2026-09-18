//! Benchmarks for `molframe-core` hot paths.
//!
//! Covers columnar chunk assembly, transparent input decompression, format
//! detection, and bond-table adjacency construction. All inputs are embedded
//! fixtures so the benches are self-contained.

use criterion::{Criterion, Throughput, black_box};
use molframe_bench::{Sample, input, input_gzip, large_cif_gz, structure};
use molframe_core::chunk::{AtomRecord, ChunkBuilder};
use molframe_core::column::Presence;
use molframe_core::element::Element;
use molframe_core::index::ResidueIndex;
use molframe_core::io::{Format, InputBuffer, Limits};
use molframe_core::optional::OptionalSymbol;
use molframe_core::symbol::{AltId, SymbolId};

/// Builds a realistic number of synthetic atom records without depending on a
/// particular structure's private fields.
fn make_records(count: usize) -> Vec<AtomRecord> {
    let mut records = Vec::with_capacity(count);
    let atom_names: [SymbolId; 4] = [
        SymbolId::from_raw(1),
        SymbolId::from_raw(2),
        SymbolId::from_raw(3),
        SymbolId::from_raw(4),
    ];
    let elements: [Element; 4] = [
        Element::NITROGEN,
        Element::CARBON,
        Element::CARBON,
        Element::OXYGEN,
    ];

    for index in 0..count {
        let slot = index % 4;
        records.push(AtomRecord {
            position: Some([as_f32(index), as_f32(index % 7), as_f32(index % 13)]),
            element: elements[slot],
            atom_name: atom_names[slot],
            auth_atom_name: OptionalSymbol::NONE,
            alternate_component_id: OptionalSymbol::NONE,
            alt_id: AltId::BLANK,
            residue: ResidueIndex::new(as_u32(index / 4)),
            occupancy: (1.0, Presence::Present),
            b_factor: (20.0, Presence::Present),
            formal_charge: (0, Presence::Inapplicable),
            atom_site_id: match as_u32(index).checked_add(1) {
                Some(value) => value,
                None => panic!("synthetic atom-site id overflowed"),
            },
        });
    }
    records
}

fn as_u32(value: usize) -> u32 {
    match u32::try_from(value) {
        Ok(value) => value,
        Err(error) => panic!("synthetic benchmark value exceeds u32: {error}"),
    }
}

fn as_f32(value: usize) -> f32 {
    let Ok(value) = u16::try_from(value) else {
        panic!("synthetic benchmark coordinate exceeds f32 fixture range")
    };
    f32::from(value)
}

fn bench_chunk_assembly(c: &mut Criterion) {
    let records = make_records(5_000);

    let mut group = c.benchmark_group("core_chunk_builder");
    group.throughput(Throughput::Elements(records.len() as u64));
    group.bench_function("push_finish", |b| {
        b.iter(|| {
            let mut builder = ChunkBuilder::new();
            builder.start_model(0);
            for record in &records {
                builder.push(*record);
            }
            black_box(builder.finish());
        });
    });
    group.finish();
}

fn bench_input_open(c: &mut Criterion) {
    let Some(medium_cif) = Sample::Medium.cif() else {
        panic!("medium sample has no CIF fixture")
    };
    let uncompressed = input(medium_cif);
    let compressed = input_gzip(large_cif_gz());

    let mut group = c.benchmark_group("core_input_open");
    group.throughput(Throughput::Bytes(uncompressed.len() as u64));
    group.bench_function("uncompressed", |b| {
        b.iter(|| {
            let _ = black_box(InputBuffer::open(
                "crates/molframe-bench/data/4hhb.cif",
                Limits::default(),
            ));
        });
    });
    group.throughput(Throughput::Bytes(large_cif_gz().len() as u64));
    group.bench_function("gzip", |b| {
        b.iter(|| {
            let _ = black_box(InputBuffer::from_reader(
                std::io::Cursor::new(large_cif_gz().to_vec()),
                Limits::default(),
            ));
        });
    });
    black_box((uncompressed.len(), compressed.len()));
    group.finish();
}

fn bench_format_detect(c: &mut Criterion) {
    let Some(small_cif) = Sample::Small.cif() else {
        panic!("small sample has no CIF fixture")
    };
    let Some(small_pdb) = Sample::Small.pdb() else {
        panic!("small sample has no PDB fixture")
    };
    let cif_input = input(small_cif);
    let pdb_input = input(small_pdb);
    let bcif_input = input(Sample::Small.bcif());

    let mut group = c.benchmark_group("core_format_detect");
    group.bench_function("cif", |b| {
        b.iter(|| {
            let _ = black_box(Format::detect(Format::Auto, &cif_input, Some("1ubq.cif")));
        });
    });
    group.bench_function("pdb", |b| {
        b.iter(|| {
            let _ = black_box(Format::detect(Format::Auto, &pdb_input, Some("1ubq.pdb")));
        });
    });
    group.bench_function("bcif", |b| {
        b.iter(|| {
            let _ = black_box(Format::detect(Format::Auto, &bcif_input, Some("1ubq.bcif")));
        });
    });
    group.finish();
}

fn bench_bond_adjacency(c: &mut Criterion) {
    for sample in Sample::MODELS {
        let structure = structure(sample);
        let bonds = structure.data().bonds.clone();
        let atom_count = structure.atom_count();
        c.bench_function(&format!("core_bond_adjacency/{}", sample.label()), |b| {
            b.iter(|| {
                black_box(bonds.adjacency(atom_count));
            });
        });
    }
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_chunk_assembly(&mut criterion);
    bench_input_open(&mut criterion);
    bench_format_detect(&mut criterion);
    bench_bond_adjacency(&mut criterion);
    criterion.final_summary();
}
