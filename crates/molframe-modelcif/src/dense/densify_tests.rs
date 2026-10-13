use crate::{ModelCif, PaeKind, lower};
use molframe_core::io::InputBuffer;
use std::fmt::Write as _;

const DEFS: &str = "data_m\nloop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n_ma_qa_metric.software_group_id\n2 PAE PAE local-pairwise 1\n#\n";

const HEAD: &str = "data_m\nloop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n_ma_qa_metric.software_group_id\n2 PAE PAE local-pairwise 1\n#\nloop_\n_ma_qa_metric_local_pairwise.ordinal_id\n_ma_qa_metric_local_pairwise.model_id\n_ma_qa_metric_local_pairwise.label_asym_id_1\n_ma_qa_metric_local_pairwise.label_seq_id_1\n_ma_qa_metric_local_pairwise.label_comp_id_1\n_ma_qa_metric_local_pairwise.label_asym_id_2\n_ma_qa_metric_local_pairwise.label_seq_id_2\n_ma_qa_metric_local_pairwise.label_comp_id_2\n_ma_qa_metric_local_pairwise.metric_id\n_ma_qa_metric_local_pairwise.metric_value\n";

fn table(values: &[&str]) -> String {
    let names = ["MET", "ALA", "GLY"];
    let size = 3;
    let mut text = String::from(HEAD);
    for (index, value) in values.iter().enumerate() {
        let (a, b) = (index / size, index % size);
        let _ = writeln!(
            text,
            "{} 1 A {} {} A {} {} 2 {value}",
            index + 1,
            a + 1,
            names[a],
            b + 1,
            names[b]
        );
    }
    text.push_str("#\n");
    text
}

fn model(text: &str) -> ModelCif {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (document, _) = molframe_cif::parse(&input).expect("parse");
    lower(&document).expect("lower").0
}

const GRID: [&str; 9] = [
    "0.25", "1.5", "31.75", "2.5", "0.25", "4.25", "9.5", "3.75", "0.5",
];

#[test]
fn exact_table_becomes_a_matrix_in_the_narrowest_exact_layout() {
    let model = model(&table(&GRID));
    let dense = model.dense_pae().expect("dense");
    assert!(model.category("ma_qa_metric_local_pairwise").is_none());
    assert_eq!(dense.matrix().size(), 3);
    assert_eq!(dense.matrix().kind(), PaeKind::U16);
    assert_eq!(dense.matrix().get_f64(0, 2), Some(31.75));
    assert_eq!(dense.metric_id(), "2");
    assert_eq!(model.confidence().pae().count(), 9);
    let last = model.confidence().pae().last().expect("pair");
    assert_eq!((last.second_sequence_id, last.value), (3, 0.5));
}

#[test]
fn pairs_follow_the_source_rows() {
    let model = model(&table(&GRID));
    let values: Vec<f64> = model.confidence().pae().map(|pair| pair.value).collect();
    let expected: Vec<f64> = GRID
        .iter()
        .map(|text| text.parse().expect("number"))
        .collect();
    assert_eq!(values, expected);
}

#[test]
fn written_matrix_reads_back_to_the_same_matrix() {
    let source = model(&table(&GRID));
    let mut text = String::new();
    crate::write::dense_pae::write_dense_pae(&mut text, source.dense_pae().expect("dense"));
    let reread = model(&format!("{DEFS}{text}"));
    assert_eq!(source.dense_pae(), reread.dense_pae());
}

#[test]
fn gaps_extra_items_and_inexact_values_stay_compact() {
    let mut gap = table(&GRID);
    gap = gap.replace("5 1 A 2 ALA A 2 ALA 2 0.25\n", "");
    assert!(model(&gap).dense_pae().is_none());
    let renamed = table(&GRID).replace(
        "metric_value",
        "metric_value\n_ma_qa_metric_local_pairwise.extra",
    );
    assert!(model(&renamed).dense_pae().is_none());
    let odd = table(&[
        "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9",
    ]);
    let dense = model(&odd);
    let dense = dense.dense_pae().expect("tenths are exact fixed point");
    assert_eq!(dense.matrix().get_f64(2, 2), Some(0.9));
}

#[test]
fn a_replacement_matrix_must_keep_the_side() {
    use crate::{PaeLayout, PaeMatrix};
    let layout = PaeLayout {
        kind: PaeKind::U8,
        denominator: 1,
        decimals: None,
    };
    let mut model = model(&table(&GRID));
    let wrong = PaeMatrix::zeroed(2, layout, 1 << 10).expect("matrix");
    assert!(model.replace_pae_matrix(wrong).is_err());
    let right = PaeMatrix::zeroed(3, layout, 1 << 10).expect("matrix");
    model.replace_pae_matrix(right).expect("same side");
    assert_eq!(
        model.dense_pae().map(|dense| dense.matrix().kind()),
        Some(PaeKind::U8)
    );
    let mut empty = ModelCif::default();
    let none = PaeMatrix::zeroed(1, layout, 1 << 10).expect("matrix");
    assert!(empty.replace_pae_matrix(none).is_err());
}
