//! Hostile containers fed to the incremental batch source: forged row counts,
//! absurd declared sizes in every column encoding, truncation and bit flips.
//! Each must end in an error or a bounded batch, never a panic, and the
//! source's reservations must stay inside the execution budget throughout.

use super::BcifBatchSource;
use crate::container::{EncodedBlock, EncodedCategory, EncodedColumn, EncodedFile};
use crate::{DataType, EncodedData, Encoding, encode_floats, encode_integers, encode_strings};
use molframe_core::{
    Backpressure, BatchDemand, BatchSource, ChunkId, DatasetId, ExecutionContext, InputBuffer,
    LogicalRow, MemoryBudget, ReadOptions, ScratchPolicy, StructureBatchError,
};
use proptest::prelude::*;

const BUDGET: usize = 4 << 20;
const COLUMNS: [&str; 11] = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_comp_id",
    "label_asym_id",
    "label_seq_id",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "pdbx_PDB_model_num",
];

fn context() -> ExecutionContext {
    ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(BUDGET).expect("non-zero budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("bounded context")
}

/// Pulls every batch; returns the rows read or the first error.
fn drive(bytes: Vec<u8>) -> Result<usize, StructureBatchError> {
    let context = context();
    let mut source = BcifBatchSource::new(
        InputBuffer::from_bytes(bytes),
        ReadOptions::default(),
        DatasetId::new(1),
        ChunkId::new(0),
        LogicalRow::new(0),
        61,
        &context,
    )?;
    let mut rows = 0;
    for _ in 0..64 {
        assert!(context.reserved_bytes() <= BUDGET);
        match source.next_batch(BatchDemand::new(1024, 1 << 20), &context)? {
            Backpressure::Ready(batch) => rows += batch.batch().positions().len(),
            Backpressure::Finished => return Ok(rows),
            Backpressure::Pending => {}
        }
    }
    Ok(rows)
}

fn valid_columns(rows: usize) -> Vec<EncodedColumn> {
    let ints = (1..=i64::try_from(rows).expect("small")).collect::<Vec<_>>();
    let floats = (1..=rows)
        .map(molframe_core::numeric::usize_to_f64)
        .collect::<Vec<_>>();
    let words = |word: &str| vec![word.to_owned(); rows];
    let text = |word: &str| column(encode_strings(&words(word)).expect("strings"));
    vec![
        text("ATOM"),
        column(encode_integers(&ints).expect("ints")),
        text("C"),
        text("CA"),
        text("GLY"),
        text("A"),
        column(encode_integers(&ints).expect("ints")),
        column(encode_floats(&floats).expect("floats")),
        column(encode_floats(&floats).expect("floats")),
        column(encode_floats(&floats).expect("floats")),
        column(encode_integers(&vec![1; rows]).expect("ints")),
    ]
}

fn column(data: EncodedData) -> EncodedColumn {
    EncodedColumn {
        name: String::new(),
        data,
        mask: None,
    }
}

fn file(declared_rows: usize, columns: Vec<EncodedColumn>) -> Vec<u8> {
    let columns = columns
        .into_iter()
        .zip(COLUMNS)
        .map(|(mut column, name)| {
            column.name = name.to_owned();
            column
        })
        .collect();
    let file = EncodedFile {
        version: "0.3.0".to_owned(),
        encoder: "hostile-test".to_owned(),
        data_blocks: vec![EncodedBlock {
            header: "test".to_owned(),
            categories: vec![EncodedCategory {
                name: "_atom_site".to_owned(),
                row_count: declared_rows,
                columns,
            }],
        }],
    };
    crate::messagepack::consume_to_vec_named(file).expect("serialize fixture")
}

fn rejected(result: Result<usize, StructureBatchError>, what: &str) {
    match result {
        Err(StructureBatchError::Diagnostic(_) | StructureBatchError::Memory(_)) => {}
        other => panic!("{what}: expected a diagnostic, got {other:?}"),
    }
}

fn hostile_encodings() -> Vec<(&'static str, EncodedData)> {
    let pairs = [7_i32.to_le_bytes(), i32::MAX.to_le_bytes()].concat();
    let run = |src_size| EncodedData {
        encoding: vec![
            Encoding::RunLength {
                src_type: DataType::Int32,
                src_size,
            },
            Encoding::ByteArray {
                r#type: DataType::Int32,
            },
        ],
        data: pairs.clone(),
    };
    let packing = |src_size| EncodedData {
        encoding: vec![
            Encoding::IntegerPacking {
                byte_count: 1,
                is_unsigned: true,
                src_size,
            },
            Encoding::ByteArray {
                r#type: DataType::Uint8,
            },
        ],
        data: vec![255, 1],
    };
    let strings = |src_size| EncodedData {
        encoding: vec![Encoding::StringArray {
            data_encoding: vec![
                Encoding::RunLength {
                    src_type: DataType::Int8,
                    src_size,
                },
                Encoding::ByteArray {
                    r#type: DataType::Int8,
                },
            ],
            string_data: "A".to_owned(),
            offset_encoding: vec![Encoding::ByteArray {
                r#type: DataType::Int8,
            }],
            offsets: vec![0, 1],
        }],
        data: vec![0, 120],
    };
    let wild_index = EncodedData {
        data: vec![200, 200, 200],
        ..strings(3)
    };
    vec![
        ("run length, usize::MAX", run(usize::MAX)),
        ("run length, 2^40", run(1 << 40)),
        ("run length, count beyond size", run(4)),
        ("integer packing, usize::MAX", packing(usize::MAX)),
        ("integer packing, 2^40", packing(1 << 40)),
        ("dictionary run, usize::MAX", strings(usize::MAX)),
        ("dictionary run, 2^40", strings(1 << 40)),
        ("dictionary index out of range", wild_index),
        (
            "byte array of ragged length",
            EncodedData {
                encoding: vec![Encoding::ByteArray {
                    r#type: DataType::Float64,
                }],
                data: vec![1, 2, 3],
            },
        ),
        (
            "fixed point over empty bytes",
            EncodedData {
                encoding: vec![
                    Encoding::FixedPoint {
                        factor: 0.0,
                        src_type: DataType::Float32,
                    },
                    Encoding::ByteArray {
                        r#type: DataType::Int32,
                    },
                ],
                data: Vec::new(),
            },
        ),
    ]
}

#[test]
fn a_forged_row_count_errors_without_reserving_for_it() {
    for declared in [4, 1 << 20, 1 << 27, usize::MAX] {
        let bytes = file(declared, valid_columns(3));
        rejected(drive(bytes), &format!("{declared} declared rows"));
    }
}

#[test]
fn a_row_count_below_the_data_is_bounded_by_the_declaration() {
    let rows = drive(file(2, valid_columns(3))).expect("declared prefix is readable");
    assert!(rows <= 3);
}

#[test]
fn absurd_declared_sizes_in_any_column_never_panic_or_overrun() {
    for (what, hostile) in hostile_encodings() {
        for target in 0..COLUMNS.len() {
            let mut columns = valid_columns(3);
            columns[target] = column(hostile.clone());
            // Columns the batch source does not consume may legitimately be
            // skipped, so success is allowed; a failure must be a diagnostic.
            if let Err(error) = drive(file(3, columns)) {
                assert!(
                    matches!(
                        error,
                        StructureBatchError::Diagnostic(_) | StructureBatchError::Memory(_)
                    ),
                    "{what} in {}: {error:?}",
                    COLUMNS[target]
                );
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(600))]

    #[test]
    fn truncation_never_panics_and_stays_in_budget(at in any::<prop::sample::Index>()) {
        let bytes = file(3, valid_columns(3));
        let cut = at.index(bytes.len() + 1);
        drop(drive(bytes[..cut].to_vec()));
    }

    #[test]
    fn byte_flips_never_panic_and_stay_in_budget(
        edits in proptest::collection::vec((any::<prop::sample::Index>(), any::<u8>()), 1..16),
    ) {
        let mut bytes = file(3, valid_columns(3));
        for (at, value) in edits {
            let at = at.index(bytes.len());
            bytes[at] = value;
        }
        drop(drive(bytes));
    }

    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..512)) {
        drop(drive(bytes));
    }
}
