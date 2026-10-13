use crate::{ModelCif, ModelCifOptions, PaeKind, lower, read_compact, read_compact_with_options};
use molframe_core::io::InputBuffer;
use std::fmt::Write as _;

const DEFS: &str = "data_m\nloop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n_ma_qa_metric.software_group_id\n2 PAE PAE local-pairwise 1\n#\n";
const ITEMS: [&str; 10] = [
    "ordinal_id",
    "model_id",
    "label_asym_id_1",
    "label_seq_id_1",
    "label_comp_id_1",
    "label_asym_id_2",
    "label_seq_id_2",
    "label_comp_id_2",
    "metric_id",
    "metric_value",
];

fn header(items: &[&str]) -> String {
    let mut text = String::from(DEFS);
    text.push_str("loop_\n");
    for item in items {
        let _ = writeln!(text, "_ma_qa_metric_local_pairwise.{item}");
    }
    text
}

fn rows(size: usize, value: impl Fn(usize, usize) -> String) -> String {
    let mut text = header(&ITEMS);
    for row in 0..size {
        for column in 0..size {
            let _ = writeln!(
                text,
                "{} 1 A {} MET A {} MET 2 {}",
                row * size + column + 1,
                row + 1,
                column + 1,
                value(row, column)
            );
        }
    }
    text.push_str("#\n");
    text
}

fn streamed(text: &str) -> ModelCif {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    read_compact(&input).expect("read").0
}

fn reference(text: &str) -> ModelCif {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (document, _) = molframe_cif::parse(&input).expect("parse");
    lower(&document).expect("lower").0
}

fn cells(model: &ModelCif) -> Vec<Vec<String>> {
    model
        .categories
        .iter()
        .map(|category| {
            (0..category.row_count())
                .flat_map(|row| {
                    category
                        .items()
                        .iter()
                        .map(move |item| (row, item.clone()))
                        .collect::<Vec<_>>()
                })
                .map(|(row, item)| {
                    let value = category.column(&item).and_then(|column| column.value(row));
                    format!("{}.{item}[{row}]={value:?}", category.name())
                })
                .collect()
        })
        .collect()
}

/// The streamed read must hold exactly what the document route holds.
fn agrees(text: &str) -> ModelCif {
    let model = streamed(text);
    let expected = reference(text);
    assert_eq!(cells(&model), cells(&expected));
    assert_eq!(model.dense_pae(), expected.dense_pae());
    model
}

fn with_extra_item(text: &str) -> String {
    let mut out = String::new();
    for line in text.lines() {
        out.push_str(line);
        if line.starts_with(|c: char| c.is_ascii_digit()) && line.matches(' ').count() > 5 {
            out.push_str(" x");
        }
        if line.ends_with("metric_value") {
            out.push_str("\n_ma_qa_metric_local_pairwise.extra");
        }
        out.push('\n');
    }
    out
}

#[test]
fn lossless_table_streams_into_a_matrix() {
    let text = rows(5, |r, c| format!("{}.5", r + c));
    let model = agrees(&text);
    let dense = model.dense_pae().expect("dense");
    assert_eq!(dense.matrix().size(), 5);
    assert_eq!(dense.matrix().get_f64(2, 4), Some(6.5));
    assert!(model.category("ma_qa_metric_local_pairwise").is_none());
}

#[test]
fn later_values_widen_the_cell_layout_without_loss() {
    let text = rows(4, |r, c| {
        if (r, c) == (3, 3) {
            "0.1234".into()
        } else {
            format!("{}", r + c)
        }
    });
    let model = agrees(&text);
    let matrix = model.dense_pae().expect("dense").matrix();
    assert_ne!(matrix.kind(), PaeKind::U8);
    assert_eq!(matrix.get_f64(3, 3), Some(0.1234));
    assert_eq!(matrix.get_f64(1, 2), Some(3.0));
}

#[test]
fn single_residue_table_is_a_one_cell_matrix() {
    let model = agrees(&rows(1, |_, _| "1.5".into()));
    assert_eq!(model.dense_pae().expect("dense").matrix().size(), 1);
}

#[test]
fn tables_that_are_not_a_lossless_fit_stay_compact_and_unchanged() {
    let good = rows(3, |r, c| format!("{}.25", r * 3 + c));
    let cases = [
        good.replace("5 1 A 2 MET A 2 MET 2 4.25\n", ""),
        good.replace("1 1 A 1 MET A 1 MET 2", "1 3 A 1 MET A 1 MET 2"),
        with_extra_item(&good),
        good.replace("A 2 MET A 2 MET 2 4.25", "A 2 MET A 2 MET 2 ?"),
        good.replace("A 3 MET A 1 MET 2 6.25", "A 3 MET A 2 MET 2 6.25"),
        good.replace("2 6.25", "2 6.2500000000001"),
        good.replace(" 2 0.25\n", " 2 .\n"),
    ];
    for text in cases {
        let model = agrees(&text);
        assert!(model.dense_pae().is_none());
        assert!(model.category("ma_qa_metric_local_pairwise").is_some());
    }
}

#[test]
fn a_metric_that_is_not_pae_keeps_the_compact_table() {
    let text = rows(3, |r, c| format!("{}", r + c)).replace("PAE PAE", "PLDDT pLDDT");
    let model = agrees(&text);
    assert!(model.dense_pae().is_none());
    assert_eq!(
        model
            .category("ma_qa_metric_local_pairwise")
            .map(crate::ModelCategory::row_count),
        Some(9)
    );
}

#[test]
fn fallback_keeps_category_position_and_row_order() {
    let table = rows(2, |r, c| format!("{}", r + c));
    let (defs, body) = table.split_at(DEFS.len());
    let text = format!(
        "{defs}{body}loop_\n_ma_software_group.ordinal_id\n_ma_software_group.group_id\n1 1\n#\n"
    )
    .replace("2 MET A 2 MET 2 2", "2 MET A 2 MET 3 2");
    let model = agrees(&text);
    let names: Vec<&str> = model
        .categories
        .iter()
        .map(crate::ModelCategory::name)
        .collect();
    assert_eq!(
        names,
        [
            "ma_qa_metric",
            "ma_qa_metric_local_pairwise",
            "ma_software_group"
        ]
    );
}

#[test]
fn dense_table_keeps_its_slot_among_other_categories() {
    let table = rows(2, |r, c| format!("{}", r + c));
    let text = format!(
        "{table}loop_\n_ma_software_group.ordinal_id\n_ma_software_group.group_id\n1 1\n#\n"
    );
    let model = agrees(&text);
    assert_eq!(model.dense_pae().map(|dense| dense.slot), Some(1));
}

#[test]
fn streaming_never_holds_the_pair_per_row_table() {
    let text = rows(120, |r, c| format!("{}", (r * 7 + c) % 31));
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    assert!(streamed(&text).dense_pae().is_some());
    let tight = ModelCifOptions::new().with_memory_limit(64 * 1024);
    let (model, _) = read_compact_with_options(&input, tight).expect("dense within a small limit");
    assert_eq!(
        model.dense_pae().map(|dense| dense.matrix().size()),
        Some(120)
    );
}

#[test]
fn widening_mid_stream_also_stays_within_a_small_limit() {
    let text = rows(120, |r, c| {
        if (r, c) == (119, 119) {
            "0.125".into()
        } else {
            format!("{}", (r + c) % 31)
        }
    });
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let tight = ModelCifOptions::new().with_memory_limit(64 * 1024);
    let (model, _) = read_compact_with_options(&input, tight).expect("dense within a small limit");
    let matrix = model.dense_pae().expect("dense").matrix();
    assert_eq!(matrix.get_f64(119, 119), Some(0.125));
    assert_eq!(matrix.get_f64(0, 5), Some(5.0));
}
