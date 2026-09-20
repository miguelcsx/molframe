//! Criterion coverage for the all-features facade read/write path.

use criterion::{Criterion, Throughput, black_box};
use molframe::{
    Cost, OperationMetadata, ReadOptions, Workflow, WorkflowInputs, read_bytes, write_mmcif,
};
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
    let structure: molframe::Structure = structure.into();
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

    let coordinates = structure.coordinates().to_vec();
    let build = || {
        let mut workflow = Workflow::new();
        let input = match workflow.input::<Vec<[f32; 3]>>("coordinates", Cost::Borrow) {
            Ok(input) => input,
            Err(error) => panic!("workflow input failed: {error}"),
        };
        let centroid = match workflow.map(
            input,
            OperationMetadata::new("geometry.centroid", Cost::Materialize)
                .with_cse_key("geometry.centroid"),
            |positions, _| {
                let mut sum = [0.0_f64; 3];
                for position in positions {
                    for axis in 0..3 {
                        sum[axis] += f64::from(position[axis]);
                    }
                }
                let Ok(count) = u32::try_from(positions.len()) else {
                    return Err(molframe::WorkflowError::Operation {
                        operation: "geometry.centroid",
                        message: "coordinate count exceeds u32".into(),
                    });
                };
                let scale = 1.0 / f64::from(count);
                Ok(sum.map(|value| value * scale))
            },
        ) {
            Ok(centroid) => centroid,
            Err(error) => panic!("workflow operation failed: {error}"),
        };
        if let Err(error) = workflow.output("centroid", centroid) {
            panic!("workflow output failed: {error}");
        }
        workflow
    };
    c.bench_function("workflow_compile", |b| {
        b.iter(|| {
            let workflow = build();
            black_box(workflow.compile())
        });
    });
    let compiled = match build().compile() {
        Ok(compiled) => compiled,
        Err(error) => panic!("workflow compilation failed: {error}"),
    };
    let mut inputs = WorkflowInputs::new();
    inputs.insert("coordinates", coordinates);
    let context = molframe::ExecutionContext::default();
    c.bench_function("workflow_reuse", |b| {
        b.iter(|| black_box(compiled.run(&inputs, &context)));
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_facade(&mut criterion);
    criterion.final_summary();
}
