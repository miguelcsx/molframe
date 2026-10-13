use super::*;
use crate::read_compact;
use molframe_cif::CifValue;
use molframe_core::io::InputBuffer;
use std::fmt::Write as _;

fn model(size: usize) -> ModelCif {
    let mut text = String::from(
        "data_m\nloop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n_ma_qa_metric.software_group_id\n2 PAE PAE local-pairwise 1\n#\nloop_\n_ma_qa_metric_local_pairwise.ordinal_id\n_ma_qa_metric_local_pairwise.model_id\n_ma_qa_metric_local_pairwise.label_asym_id_1\n_ma_qa_metric_local_pairwise.label_seq_id_1\n_ma_qa_metric_local_pairwise.label_comp_id_1\n_ma_qa_metric_local_pairwise.label_asym_id_2\n_ma_qa_metric_local_pairwise.label_seq_id_2\n_ma_qa_metric_local_pairwise.label_comp_id_2\n_ma_qa_metric_local_pairwise.metric_id\n_ma_qa_metric_local_pairwise.metric_value\n",
    );
    for row in 0..size {
        for column in 0..size {
            let _ = writeln!(
                text,
                "{} 1 A {} MET B {} GLY 2 {}.5",
                row * size + column + 1,
                row + 1,
                column + 1,
                (row + column) % 30
            );
        }
    }
    text.push_str("#\n");
    let input = InputBuffer::from_bytes(text.into_bytes());
    read_compact(&input).expect("read").0
}

#[test]
fn export_cells_agree_with_the_rendered_text() {
    let model = model(5);
    let mut text = b"data_x\n".to_vec();
    crate::write_extension_to(&model, &mut text).expect("renders");
    let (document, _) = molframe_cif::parse(&InputBuffer::from_bytes(text)).expect("parses");
    let tables = model.export_tables();
    assert!(!tables.is_empty());
    let mut checked = 0;
    for table in &tables {
        let category = document
            .blocks()
            .flat_map(molframe_cif::DataBlock::categories)
            .find(|category| category.name() == table.name())
            .expect("category rendered");
        assert_eq!(category.row_count(), table.row_count());
        for column in 0..table.column_count() {
            let values = category.column(table.item(column)).expect("column");
            for row in 0..table.row_count() {
                let expected = values.get(row).expect("value");
                let agrees = match (table.cell(column, row), expected) {
                    (ExportCell::Integer(a), CifValue::Integer(b)) => a == *b,
                    (ExportCell::Float(a), CifValue::Float(b)) => (a - b).abs() < f64::EPSILON,
                    (ExportCell::Text(a), CifValue::Text(b)) => a == &**b,
                    (ExportCell::Unknown, CifValue::Unknown)
                    | (ExportCell::Inapplicable, CifValue::Inapplicable) => true,
                    _ => false,
                };
                assert!(agrees, "{} {column} {row}", table.name());
                checked += 1;
            }
        }
    }
    assert!(checked > 5 * 5 * 10);
}

#[test]
fn dense_table_reports_one_row_per_pair_and_is_read_only_borrow() {
    let model = model(4);
    let dense = model
        .export_tables()
        .into_iter()
        .find(|table| table.name() == "ma_qa_metric_local_pairwise")
        .expect("dense table");
    assert_eq!(dense.row_count(), 16);
    assert_eq!(dense.cell(0, 15), ExportCell::Integer(16));
    assert_eq!(dense.cell(99, 0), ExportCell::Unknown);
}
