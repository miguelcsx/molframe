use super::category::ModelBuilder;
use super::{read_compact, read_compact_with_options};
use crate::lower::tests::MODEL_CIF;
use crate::{ModelCif, ModelCifError, ModelCifOptions, ModelCifReadError, lower};
use pdbiox_cif::CifScalar;
use pdbiox_core::io::InputBuffer;

#[test]
fn direct_projection_is_cell_exact_with_document_lowering() {
    let input = InputBuffer::from_bytes(MODEL_CIF.as_bytes().to_vec());
    let (direct, direct_findings) = read_compact(&input).expect("direct projection succeeds");
    let (document, parse_findings) = pdbiox_cif::parse(&input).expect("document parse succeeds");
    let (reference, reference_findings) = lower(&document).expect("document lowering succeeds");
    assert_eq!(direct_findings, parse_findings);
    assert_eq!(direct_findings, reference_findings);
    assert_cells_equal(&direct, &reference);
}

#[test]
fn adaptive_columns_preserve_types_sentinels_and_widening_exactly() {
    let source = "data_model\nloop_\n\
_ma_future.id\n_ma_future.integer\n_ma_future.numeric\n_ma_future.mixed\n\
1 1 1 1\n\
2 200 200 200\n\
3 40000 40000 40000\n\
4 5000000000 5000000000 5000000000\n\
5 . 1.5 1.5\n\
6 ? ? text\n\
7 -7 . .\n\
8 8 ? ?\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (direct, direct_findings) = read_compact(&input).expect("direct projection succeeds");
    let (document, parse_findings) = pdbiox_cif::parse(&input).expect("document parse succeeds");
    let (reference, reference_findings) = lower(&document).expect("document lowering succeeds");
    assert_eq!(direct_findings, parse_findings);
    assert_eq!(direct_findings, reference_findings);
    assert_cells_equal(&direct, &reference);
}

#[test]
fn direct_and_document_paths_report_identical_semantic_findings() {
    let source = "data_model\nloop_\n\
_ma_qa_metric_local.model_id\n_ma_qa_metric_local.label_asym_id\n\
_ma_qa_metric_local.label_seq_id\n_ma_qa_metric_local.metric_id\n\
_ma_qa_metric_local.metric_value\n1 A nope 2 infinite\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (direct, direct_findings) = read_compact(&input).expect("direct projection succeeds");
    let (document, parse_findings) = pdbiox_cif::parse(&input).expect("document parse succeeds");
    let (reference, reference_findings) = lower(&document).expect("document lowering succeeds");
    assert_cells_equal(&direct, &reference);
    assert_eq!(parse_findings.len(), 0);
    assert_eq!(direct_findings, reference_findings);
}

#[test]
fn direct_projection_enforces_the_same_memory_policy() {
    let input = InputBuffer::from_bytes(MODEL_CIF.as_bytes().to_vec());
    let error = read_compact_with_options(&input, ModelCifOptions::new().with_memory_limit(1))
        .expect_err("one byte cannot retain the fixture");
    assert!(matches!(
        error,
        ModelCifReadError::Projection(crate::ModelCifError::MemoryLimit { limit: 1, .. })
    ));
}

#[test]
fn memory_policy_counts_the_temporary_integer_widening_buffer() {
    // A budget large enough that the measurement itself is never the bound.
    let mut measured = integer_builder(500_000_000);
    measured
        .push("ma_future", "value", CifScalar::Integer(128))
        .expect("measurement builder widens");
    let retained_after_widening = measured.tracked_bytes();

    let mut bounded = integer_builder(retained_after_widening);
    let error = bounded
        .push("ma_future", "value", CifScalar::Integer(128))
        .expect_err("the transient old and widened buffers exceed retained bytes");
    let ModelCifError::MemoryLimit { required, limit } = error else {
        panic!("expected the peak-memory policy")
    };
    assert_eq!(limit, retained_after_widening);
    assert!(required > retained_after_widening);
}

#[test]
fn syntax_failure_is_not_relabelled_as_a_projection_failure() {
    let input = InputBuffer::from_bytes(Vec::new());
    let error = read_compact(&input).expect_err("an empty document is invalid");
    let ModelCifReadError::Syntax(findings) = error else {
        panic!("expected syntax findings");
    };
    assert_eq!(findings.len(), 1);
    assert_eq!(findings[0].code(), pdbiox_core::Code::E1106);
}

fn assert_cells_equal(actual: &ModelCif, expected: &ModelCif) {
    assert_eq!(actual.categories.len(), expected.categories.len());
    for (actual, expected) in actual.categories.iter().zip(&expected.categories) {
        assert_eq!(actual.name(), expected.name());
        assert_eq!(actual.items(), expected.items());
        assert_eq!(actual.row_count(), expected.row_count());
        for item in expected.items() {
            for row in 0..expected.row_count() {
                assert_eq!(actual.value(item, row), expected.value(item, row));
            }
        }
    }
}

fn integer_builder(limit: usize) -> ModelBuilder {
    let mut builder = ModelBuilder::new(limit);
    for value in 0_i64..64 {
        builder
            .push("ma_future", "value", CifScalar::Integer(value))
            .expect("fixture stays inside its limit");
    }
    builder
}
