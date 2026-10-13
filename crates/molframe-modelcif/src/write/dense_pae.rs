//! Streams a dense PAE matrix back as the pair-per-row `ModelCIF` table.
//!
//! Rows are produced one at a time from the matrix, so writing never
//! materialises the table the matrix replaced.

use crate::DensePae;
use molframe_cif::write_quoted;
use std::fmt;

pub(crate) fn write_dense_pae(output: &mut impl fmt::Write, dense: &DensePae) {
    let size = dense.matrix().size();
    if size == 0 {
        return;
    }
    let _ = output.write_str("loop_\n");
    for item in dense.items() {
        let _ = writeln!(output, "_ma_qa_metric_local_pairwise.{item}");
    }
    for row in 0..size {
        for column in 0..size {
            write_row(output, dense, row, column);
        }
    }
    let _ = output.write_str("#\n");
}

fn write_row(output: &mut impl fmt::Write, dense: &DensePae, row: usize, column: usize) {
    let size = dense.matrix().size();
    let (Some(first), Some(second)) = (dense.residues().get(row), dense.residues().get(column))
    else {
        return;
    };
    for (index, item) in dense.items().iter().enumerate() {
        if index > 0 {
            let _ = output.write_char(' ');
        }
        match &**item {
            "ordinal_id" => {
                let _ = write!(output, "{}", row * size + column + 1);
            }
            "model_id" => write_quoted(output, dense.model_id()),
            "metric_id" => write_quoted(output, dense.metric_id()),
            "label_asym_id_1" => write_quoted(output, first.chain),
            "label_asym_id_2" => write_quoted(output, second.chain),
            "label_seq_id_1" => {
                let _ = write!(output, "{}", first.sequence);
            }
            "label_seq_id_2" => {
                let _ = write!(output, "{}", second.sequence);
            }
            "label_comp_id_1" => write_component(output, first.component),
            "label_comp_id_2" => write_component(output, second.component),
            _ => match dense.matrix().get_f64(row, column) {
                Some(value) => {
                    let _ = write!(output, "{value}");
                }
                None => {
                    let _ = output.write_char('?');
                }
            },
        }
    }
    let _ = output.write_char('\n');
}

fn write_component(output: &mut impl fmt::Write, component: Option<&str>) {
    match component {
        Some(name) => write_quoted(output, name),
        None => {
            let _ = output.write_char('?');
        }
    }
}
