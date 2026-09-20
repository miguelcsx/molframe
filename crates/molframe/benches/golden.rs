//! Criterion measurements for the executable golden workflows.

use criterion::{Criterion, Throughput, black_box};
use molframe::formats::cif::CifWriteOptions;
use molframe::formats::pdb::PdbOptions;
use molframe::geometry::{Rigid, superpose};
use molframe::{
    AltlocPolicy, AnalysisPolicy, Code, ReadOptions, read_bytes, transform, write_bcif,
    write_mmcif_with_options, write_pdb,
};
use molframe_bench::{Sample, structure_from_cif};
use molframe_core::selection::AtomSelection;

#[path = "golden/native.rs"]
mod native;
#[path = "golden/native_pending.rs"]
mod native_pending;

const MULTI_MODEL_PDB: &str = "\
HEADER    GOLDEN                              01-JAN-00   1ABC
MODEL        7
ATOM      1  N   GLY A  10A     10.000  11.000  12.000  1.00 10.00           N
ATOM      2  CA  GLY A  10A     11.000  11.000  12.000  1.00 11.00           C
ENDMDL
MODEL       11
ATOM      1  N   GLY A  10A     12.000  13.000  14.000  1.00 10.00           N
ATOM      2  CA  GLY A  10A     13.000  13.000  14.000  1.00 11.00           C
ENDMDL
END
";

const NAMESPACED_CIF: &str = "\
data_golden
loop_
_entity.id
_entity.type
1 polymer
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY LONG 1 1 0 0 0 1 10 42 GLY A N 1
ATOM 2 C CA A GLY LONG 1 1 1 0 0 0.4 11 42 GLY A CA 1
ATOM 3 C CA B GLY LONG 1 1 2 0 0 0.6 12 42 GLY A CA 1
";

fn fixture(text: &str, name: &str) -> molframe::Structure {
    match read_bytes(text.as_bytes().to_vec(), Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("golden benchmark fixture failed: {findings:?}"),
    }
}

fn medium_cif() -> (&'static [u8], molframe::Structure) {
    let Some(bytes) = Sample::Medium.cif() else {
        panic!("medium CIF fixture is required")
    };
    let Some(structure) = structure_from_cif(Sample::Medium) else {
        panic!("medium CIF structure is required")
    };
    (bytes, structure.into())
}

fn bench_gw_001(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let (bytes, _) = medium_cif();
    group.throughput(Throughput::Bytes(bytes.len() as u64));
    group.bench_function("GW-001", |b| {
        b.iter(|| {
            let result = read_bytes(bytes.to_vec(), Some("4hhb.cif"), &ReadOptions::new());
            let Ok((structure, _)) = result else {
                panic!("GW-001 read failed")
            };
            black_box((
                structure.model_count(),
                structure.chain_count(),
                structure.engine().entity_count(),
                structure.atom_count(),
            ));
        });
    });
}

fn bench_gw_002(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let bytes = NAMESPACED_CIF.as_bytes();
    group.throughput(Throughput::Bytes(bytes.len() as u64));
    group.bench_function("GW-002", |b| {
        b.iter(|| {
            let Ok((source, _)) =
                read_bytes(bytes.to_vec(), Some("golden.cif"), &ReadOptions::new())
            else {
                panic!("GW-002 source read failed")
            };
            let encoded = match write_bcif(&source) {
                Ok(encoded) => encoded,
                Err(findings) => panic!("GW-002 BinaryCIF write failed: {findings:?}"),
            };
            let Ok((round_tripped, _)) =
                read_bytes(encoded, Some("golden.bcif"), &ReadOptions::new())
            else {
                panic!("GW-002 BinaryCIF read failed")
            };
            black_box((
                round_tripped.atom_count(),
                round_tripped.coordinates().len(),
            ));
        });
    });
}

fn bench_gw_003(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    group.throughput(Throughput::Bytes(MULTI_MODEL_PDB.len() as u64));
    group.bench_function("GW-003", |b| {
        b.iter(|| {
            let Ok((structure, _)) = read_bytes(
                MULTI_MODEL_PDB.as_bytes().to_vec(),
                Some("golden.pdb"),
                &ReadOptions::new(),
            ) else {
                panic!("GW-003 PDB read failed")
            };
            let options = CifWriteOptions::new().with_block_id("1ABC");
            let rendered = match write_mmcif_with_options(&structure, &options) {
                Ok(rendered) => rendered,
                Err(error) => panic!("GW-003 mmCIF write failed: {error}"),
            };
            let round_tripped = match read_bytes(
                rendered.into_bytes(),
                Some("golden.cif"),
                &ReadOptions::new(),
            ) {
                Ok((structure, _)) => structure,
                Err(findings) => panic!("GW-003 round-trip failed: {findings:?}"),
            };
            black_box(round_tripped.atom_count());
        });
    });
}

fn bench_gw_005(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let source = NAMESPACED_CIF.replace(" 42 GLY A ", " 42 GLY LONG ");
    let structure = fixture(&source, "golden.cif");
    let result = write_pdb(&structure, &PdbOptions::new());
    assert_eq!(
        result
            .err()
            .and_then(|findings| findings.first().map(molframe::Diagnostic::code)),
        Some(Code::E4102)
    );
    group.throughput(Throughput::Bytes(source.len() as u64));
    group.bench_function("GW-005", |b| {
        b.iter(|| {
            let structure = fixture(&source, "golden.cif");
            let result = write_pdb(&structure, &PdbOptions::new());
            let _ = black_box(result);
        });
    });
}

fn bench_gw_007(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = fixture(NAMESPACED_CIF, "golden.cif");
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-007", |b| {
        b.iter(|| {
            let counts: Vec<_> = [
                AltlocPolicy::KeepAll,
                AltlocPolicy::First,
                AltlocPolicy::HighestOccupancyPerResidue,
                AltlocPolicy::HighestOccupancyPerAtom,
                AltlocPolicy::ConformerConsistent,
            ]
            .into_iter()
            .map(|altloc| {
                structure
                    .engine()
                    .resolve_altlocs(&AnalysisPolicy::default().with_altloc(altloc))
                    .value
                    .len()
            })
            .collect();
            black_box(counts);
        });
    });
}

fn bench_gw_009(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = fixture(NAMESPACED_CIF, "golden.cif");
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-009", |b| {
        b.iter(|| {
            let chain = structure.chain_at(0);
            black_box((
                chain.and_then(molframe::ChainRef::label),
                chain.and_then(molframe::ChainRef::auth_label),
            ));
        });
    });
}

fn bench_gw_011(group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let fixed = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mobile = [[2.0, 3.0, 4.0], [3.0, 3.0, 4.0], [2.0, 4.0, 4.0]];
    let fit = match superpose(&mobile, &fixed) {
        Ok(fit) => fit,
        Err(error) => panic!("GW-011 setup failed: {error:?}"),
    };
    assert!(fit.rmsd < 1.0e-6);
    let structure = fixture(NAMESPACED_CIF, "golden.cif");
    let selection = AtomSelection::All(structure.atom_count());
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-011", |b| {
        b.iter(|| {
            let fitted = match superpose(&mobile, &fixed) {
                Ok(fitted) => fitted,
                Err(error) => panic!("GW-011 fit failed: {error:?}"),
            };
            let moved =
                match transform(&structure, &selection, &Rigid::translation([1.0, 2.0, 3.0])) {
                    Ok(moved) => moved,
                    Err(findings) => panic!("GW-011 transform failed: {findings:?}"),
                };
            black_box((fitted.rmsd, moved.engine().generation().get()));
        });
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    let mut group = criterion.benchmark_group("golden_workflows");
    bench_gw_001(&mut group);
    bench_gw_002(&mut group);
    bench_gw_003(&mut group);
    bench_gw_005(&mut group);
    bench_gw_007(&mut group);
    bench_gw_009(&mut group);
    bench_gw_011(&mut group);
    native::register(&mut group);
    native_pending::register(&mut group);
    group.finish();
    criterion.final_summary();
}
