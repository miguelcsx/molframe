use super::*;
use crate::writer::document::encode_category;
use molframe_cif::{CifValue, DataBlock, Document};
use molframe_core::span::ByteSpan;

struct Grid {
    columns: Vec<(&'static str, Vec<CifValue>)>,
}

impl ColumnTable for Grid {
    fn name(&self) -> &'static str {
        "grid"
    }
    fn row_count(&self) -> usize {
        self.columns.first().map_or(0, |column| column.1.len())
    }
    fn column_count(&self) -> usize {
        self.columns.len()
    }
    fn item(&self, column: usize) -> &str {
        self.columns.get(column).map_or("", |column| column.0)
    }
    fn cell(&self, column: usize, row: usize) -> Cell<'_> {
        match self
            .columns
            .get(column)
            .and_then(|column| column.1.get(row))
        {
            Some(CifValue::Integer(value)) => Cell::Integer(*value),
            Some(CifValue::Float(value)) => Cell::Float(*value),
            Some(CifValue::Text(value)) => Cell::Text(value),
            Some(CifValue::Inapplicable) => Cell::Inapplicable,
            _ => Cell::Unknown,
        }
    }
}

fn text(value: &str) -> CifValue {
    CifValue::Text(value.into())
}

fn sample() -> Grid {
    Grid {
        columns: vec![
            (
                "id",
                vec![
                    CifValue::Integer(1),
                    CifValue::Integer(2),
                    CifValue::Integer(3),
                ],
            ),
            (
                "score",
                vec![
                    CifValue::Float(0.5),
                    CifValue::Unknown,
                    CifValue::Integer(4),
                ],
            ),
            ("chain", vec![text("A"), text("A"), CifValue::Inapplicable]),
            (
                "mixed",
                vec![CifValue::Integer(1), text("(1-60)"), CifValue::Float(2.5)],
            ),
            ("empty", vec![CifValue::Unknown; 3]),
        ],
    }
}

#[test]
fn direct_columns_match_the_document_encoder_byte_for_byte() {
    let grid = sample();
    let mut block = DataBlock::new("t");
    let category = block.category_mut("grid", ByteSpan::default());
    for (item, values) in &grid.columns {
        for value in values {
            category
                .column_mut(item)
                .push(value.clone(), molframe_cif::lexer::Quoting::Bare);
        }
    }
    let mut document = Document::new();
    document.push(block);
    let expected = document
        .blocks()
        .flat_map(molframe_cif::DataBlock::categories)
        .map(|category| encode_category(category).expect("document encodes"))
        .next()
        .expect("category");
    let actual = encode_table(&grid).expect("table encodes");
    assert_eq!(format!("{actual:?}"), format!("{expected:?}"));
}

#[test]
fn out_of_range_cells_do_not_panic() {
    let grid = sample();
    assert_eq!(grid.cell(9, 9), Cell::Unknown);
}
