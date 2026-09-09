//! Criterion coverage for real `ModelCIF` lowering and canonical writing.

use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, structure};
use pdbiox_cif::{lexer::Lexer, parse};
use pdbiox_core::io::InputBuffer;
use pdbiox_modelcif::{lower, read_compact, write_canonical};
use std::path::{Path, PathBuf};

const SMALL_MODEL_CIF: &str = "data_model\n\
loop_\n_ma_model_list.ordinal_id\n_ma_model_list.model_name\n_ma_model_list.model_type\n\
1 prediction 'Ab initio model'\n#\n\
loop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n\
1 pLDDT pLDDT local\n2 PAE PAE local-pairwise\n#\n\
loop_\n_ma_qa_metric_local.model_id\n_ma_qa_metric_local.label_asym_id\n\
_ma_qa_metric_local.label_seq_id\n_ma_qa_metric_local.metric_id\n\
_ma_qa_metric_local.metric_value\n1 A 1 1 91.5\n#\n\
loop_\n_ma_qa_metric_local_pairwise.model_id\n\
_ma_qa_metric_local_pairwise.label_asym_id_1\n\
_ma_qa_metric_local_pairwise.label_seq_id_1\n\
_ma_qa_metric_local_pairwise.label_asym_id_2\n\
_ma_qa_metric_local_pairwise.label_seq_id_2\n\
_ma_qa_metric_local_pairwise.metric_id\n\
_ma_qa_metric_local_pairwise.metric_value\n1 A 1 A 2 2 3.2\n";

fn parsed(bytes: Vec<u8>) -> pdbiox_cif::Document {
    let buffer = InputBuffer::from_bytes(bytes);
    match parse(&buffer) {
        Ok((document, _)) => document,
        Err(findings) => panic!("ModelCIF bench fixture failed: {findings:?}"),
    }
}

fn bench_lower(c: &mut Criterion) {
    let document = parsed(SMALL_MODEL_CIF.as_bytes().to_vec());
    c.bench_function("modelcif_lower/typed_small", |b| {
        b.iter(|| black_box(lower(&document)));
    });
    let compact_input = InputBuffer::from_bytes(SMALL_MODEL_CIF.as_bytes().to_vec());
    c.bench_function("modelcif_lex/typed_small", |b| {
        b.iter(|| black_box(lex_all(&compact_input)));
    });
    c.bench_function("modelcif_read_compact/typed_small", |b| {
        b.iter(|| black_box(read_compact(&compact_input)));
    });

    let path = stress_path("ma-kvko-prsa-005_qa_metrics.cif");
    if path.is_file() {
        let bytes = match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) => panic!("{} cannot be read: {error}", path.display()),
        };
        let stress = InputBuffer::from_bytes(bytes);
        let mut lex_group = c.benchmark_group("modelcif_lex_stress");
        lex_group.sample_size(10);
        lex_group.bench_function("pairwise_1946880", |b| {
            b.iter(|| black_box(lex_all(&stress)));
        });
        lex_group.finish();
        let mut group = c.benchmark_group("modelcif_read_compact_stress");
        group.sample_size(10);
        group.bench_function("pairwise_1946880", |b| {
            b.iter(|| black_box(read_compact(&stress)));
        });
        group.finish();
    }
}

fn lex_all(input: &InputBuffer) -> Result<usize, pdbiox_cif::lexer::LexError> {
    let mut lexer = Lexer::new(input.as_bytes())?;
    let mut count = 0usize;
    while lexer.next_token()?.is_some() {
        count = count.saturating_add(1);
    }
    Ok(count)
}

fn bench_write(c: &mut Criterion) {
    let structure = structure(Sample::Tiny);
    let document = parsed(SMALL_MODEL_CIF.as_bytes().to_vec());
    let model_cif = match lower(&document) {
        Ok((model, findings)) if findings.is_empty() => model,
        Ok((_, findings)) => panic!("ModelCIF write fixture is invalid: {findings:?}"),
        Err(error) => panic!("ModelCIF write fixture exceeded resources: {error}"),
    };
    c.bench_function("modelcif_write/typed_small", |b| {
        b.iter(|| black_box(write_canonical(&structure, &model_cif)));
    });
}

fn stress_path(name: &str) -> PathBuf {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .map_or_else(PathBuf::new, Path::to_path_buf);
    root.join("target/stress-data/modelcif").join(name)
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_lower(&mut criterion);
    bench_write(&mut criterion);
    criterion.final_summary();
}
