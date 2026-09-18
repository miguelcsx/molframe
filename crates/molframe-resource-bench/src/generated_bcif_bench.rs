//! Logical tera-scale `BinaryCIF` rows encoded as bounded run-length payloads.

use super::generated_structure_bench::drain;
use super::{ResourceRecord, measure_case};
use molframe::ReadOptions;
use molframe::bcif::{DataType, EncodedData, Encoding};
use molframe::core::{
    ChunkId, DatasetId, ExecutionContext, InputBuffer, LogicalRow, ScratchPolicy,
};
use serde::Serialize;

const LOGICAL_ROW_BYTES: u64 = 49;

pub(super) fn run(minimum_logical_bytes: u64) -> Result<ResourceRecord, String> {
    let rows = minimum_logical_bytes
        .checked_add(LOGICAL_ROW_BYTES - 1)
        .and_then(|value| value.checked_div(LOGICAL_ROW_BYTES))
        .map_or(1, |value| value.max(1));
    let bytes = encoded_file(rows)?;
    measure_case("generated_bcif_batches", move || {
        let context = ExecutionContext::builder()
            .scratch_policy(ScratchPolicy::new(0))
            .build()
            .map_err(|error| format!("generated_bcif_batches: context failed: {error}"))?;
        let source = molframe::bcif::BcifBatchSource::new(
            InputBuffer::from_bytes(bytes),
            ReadOptions::new(),
            DatasetId::new(0),
            ChunkId::new(0),
            LogicalRow::new(0),
            64 * 1024,
            &context,
        )
        .map_err(|error| format!("generated_bcif_batches: open failed: {error}"))?;
        drain(source, &context, "generated_bcif_batches")
    })
}

#[derive(Serialize)]
struct File {
    version: &'static str,
    encoder: &'static str,
    #[serde(rename = "dataBlocks")]
    data_blocks: Vec<Block>,
}

#[derive(Serialize)]
struct Block {
    header: &'static str,
    categories: Vec<Category>,
}

#[derive(Serialize)]
struct Category {
    name: &'static str,
    #[serde(rename = "rowCount")]
    row_count: usize,
    columns: Vec<Column>,
}

#[derive(Serialize)]
struct Column {
    name: &'static str,
    data: EncodedData,
}

fn encoded_file(rows: u64) -> Result<Vec<u8>, String> {
    let row_count = usize::try_from(rows)
        .map_err(|_| "generated BinaryCIF row count exceeds usize".to_owned())?;
    let columns = vec![
        integer_column("id", 1, rows)?,
        integer_column("Cartn_x", 1, rows)?,
        integer_column("Cartn_y", 2, rows)?,
        integer_column("Cartn_z", 3, rows)?,
        integer_column("pdbx_PDB_model_num", 1, rows)?,
    ];
    let file = File {
        version: "0.3.0",
        encoder: "molframe-tera-gate",
        data_blocks: vec![Block {
            header: "tera",
            categories: vec![Category {
                name: "_atom_site",
                row_count,
                columns,
            }],
        }],
    };
    rmp_serde::to_vec_named(&file)
        .map_err(|error| format!("generated BinaryCIF serialization failed: {error}"))
}

fn integer_column(name: &'static str, value: i32, rows: u64) -> Result<Column, String> {
    let mut remaining = rows;
    let mut data = Vec::new();
    while remaining > 0 {
        let count = remaining.min(i32::MAX as u64);
        data.extend_from_slice(&value.to_le_bytes());
        data.extend_from_slice(
            &i32::try_from(count)
                .map_err(|_| "generated BinaryCIF run exceeds i32".to_owned())?
                .to_le_bytes(),
        );
        remaining -= count;
    }
    Ok(Column {
        name,
        data: EncodedData {
            encoding: vec![
                Encoding::RunLength {
                    src_type: DataType::Int32,
                    src_size: usize::try_from(rows)
                        .map_err(|_| "generated BinaryCIF source size exceeds usize".to_owned())?,
                },
                Encoding::ByteArray {
                    r#type: DataType::Int32,
                },
            ],
            data,
        },
    })
}

#[cfg(test)]
#[path = "generated_bcif_tests.rs"]
mod tests;
