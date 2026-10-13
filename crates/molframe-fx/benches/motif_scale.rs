//! Criterion coverage for motif mapping at assembly scale: plan build,
//! counting and streaming over a synthetic structure of up to 100k residues.

use criterion::{BatchSize, Criterion, black_box};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::execution::CancellationToken;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_fx::{
    AtomSite, ComponentRole, ComponentSpec, Constraint, MappingPlan, Motif, NamedConstraint,
};
use std::fmt::Write as _;
use std::sync::Arc;

/// Residues scattered uniformly in a cube of `edge` angstroms.
fn scatter(residues: u32, edge: f64) -> Structure {
    let mut source = String::from(
        "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let mut state = 12_345_u64;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f64::from(u32::try_from(state >> 40).unwrap_or(0)) / 16_777_216.0 * edge
    };
    let mut serial = 0;
    for residue in 1..=residues {
        let (name, atom, element) = match residue % 3 {
            0 => ("SER", "OG", "O"),
            1 => ("TRP", "NE1", "N"),
            _ => ("MET", "SD", "S"),
        };
        let (x, y, z) = (next(), next(), next());
        serial += 1;
        let _ = writeln!(
            source,
            "ATOM {serial} {element} {atom} {name} A {residue} {x:.3} {y:.3} {z:.3}"
        );
    }
    molframe_cif::read(
        &InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .map_or_else(
        |findings| panic!("fixture failed: {findings:?}"),
        |(structure, _)| structure,
    )
}

fn chain_motif() -> Motif {
    let parts = [("a", "SER", "OG"), ("b", "TRP", "NE1"), ("c", "MET", "SD")];
    let specs = parts.iter().map(|(name, id, atom)| {
        (
            Box::<str>::from(*name),
            ComponentSpec::new(ComponentRole::Residue)
                .component(*id)
                .require(*atom),
        )
    });
    let constraint = |name: &str, first: usize, second: usize| NamedConstraint {
        name: name.into(),
        constraint: Constraint::Distance {
            first: AtomSite::new(parts[first].0, parts[first].2),
            second: AtomSite::new(parts[second].0, parts[second].2),
            target: 4.0,
            tolerance: 1.0,
        },
    };
    Motif::new(specs, [constraint("d0", 0, 1), constraint("d1", 1, 2)])
        .unwrap_or_else(|error| panic!("motif failed: {error}"))
}

fn build(structure: &Structure, motif: &Motif) -> Arc<MappingPlan> {
    Arc::new(
        MappingPlan::build(structure, motif, None, &AnalysisPolicy::default(), true)
            .unwrap_or_else(|error| panic!("plan failed: {error}")),
    )
}

fn stream_all(plan: &Arc<MappingPlan>) -> u64 {
    let cancel = CancellationToken::new();
    let mut stream = plan.stream();
    let mut tuple = vec![0u32; plan.width()];
    let mut seen = 0u64;
    while stream
        .next_into(&mut tuple, &cancel)
        .unwrap_or_else(|error| panic!("stream failed: {error}"))
    {
        seen += 1;
    }
    seen
}

/// Fastest of several runs: robust to a loaded machine, unlike a mean.
fn report_minimum(residues: u32, structure: &Structure, motif: &Motif, plan: &Arc<MappingPlan>) {
    let cancel = CancellationToken::new();
    let fastest = |mut run: Box<dyn FnMut() + '_>| {
        (0..25)
            .map(|_| {
                let start = std::time::Instant::now();
                run();
                start.elapsed()
            })
            .min()
            .unwrap_or_default()
    };
    let built = fastest(Box::new(|| drop(black_box(build(structure, motif)))));
    let counted = fastest(Box::new(|| drop(black_box(plan.count(&cancel)))));
    let streamed = fastest(Box::new(|| {
        black_box(stream_all(plan));
    }));
    eprintln!("min {residues}: build {built:?} count {counted:?} stream {streamed:?}");
}

fn bench_scale(c: &mut Criterion) {
    let motif = chain_motif();
    let mut group = c.benchmark_group("fx_motif_scale");
    group.sample_size(10);
    for (residues, edge) in [(10_000_u32, 70.0), (100_000, 150.0)] {
        let structure = scatter(residues, edge);
        let plan = build(&structure, &motif);
        let cancel = CancellationToken::new();
        let total = plan
            .count(&cancel)
            .unwrap_or_else(|error| panic!("count failed: {error}"));
        eprintln!("residues {residues}: {total} mappings");
        report_minimum(residues, &structure, &motif, &plan);
        group.bench_function(format!("build_{residues}"), |b| {
            b.iter_batched(
                || (),
                |()| black_box(build(&structure, &motif)),
                BatchSize::PerIteration,
            );
        });
        group.bench_function(format!("count_{residues}"), |b| {
            b.iter(|| black_box(plan.count(&cancel)));
        });
        group.bench_function(format!("stream_{residues}"), |b| {
            b.iter(|| black_box(stream_all(&plan)));
        });
    }
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_scale(&mut criterion);
    criterion.final_summary();
}
