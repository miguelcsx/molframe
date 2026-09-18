use super::BcifBatchSource;
use crate::container::{EncodedBlock, EncodedCategory, EncodedColumn, EncodedFile};
use crate::{DataType, EncodedData, Encoding, encode_floats, encode_integers, encode_strings};
use molframe_core::{
    Backpressure, Batch, BatchDemand, BatchSource, ChunkId, DatasetId, ExecutionContext,
    InputBuffer, LogicalRow, MemoryBudget, ReadOptions, ScratchPolicy, StructureBatchError,
};

#[test]
fn batches_decode_string_numeric_mask_and_coordinate_columns_incrementally() {
    let bytes = atom_file();
    let context = ExecutionContext::default();
    let source = InputBuffer::from_bytes(bytes);
    let mut batches = BcifBatchSource::new(
        source,
        ReadOptions::default(),
        DatasetId::new(7),
        ChunkId::new(9),
        LogicalRow::new(11),
        37,
        &context,
    )
    .expect("batch source");
    let pending = batches
        .next_batch(BatchDemand::new(0, 0), &context)
        .expect("zero demand");
    assert!(matches!(pending, Backpressure::Pending));
    let ready = batches
        .next_batch(BatchDemand::new(2, 256 * 1024), &context)
        .expect("first batch");
    let Backpressure::Ready(first) = ready else {
        panic!("first batch was not ready");
    };
    assert_eq!(first.batch().descriptor().dataset(), DatasetId::new(7));
    assert_eq!(
        first.batch().descriptor().logical_start(),
        LogicalRow::new(11)
    );
    assert_eq!(first.batch().models(), &[1, 1]);
    assert_eq!(
        first.batch().positions(),
        &[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    );
    let atom = first.batch().dictionary().resolve(first.batch().atoms()[0]);
    assert_eq!(atom, Some("CA"));
    drop(first);
    let ready = batches
        .next_batch(BatchDemand::new(2, 256 * 1024), &context)
        .expect("second batch");
    let Backpressure::Ready(second) = ready else {
        panic!("second batch was not ready");
    };
    assert_eq!(
        second.batch().descriptor().logical_start(),
        LogicalRow::new(13)
    );
    assert_eq!(second.batch().positions(), &[[7.0, 8.0, 9.0]]);
    drop(second);
    assert!(matches!(
        batches
            .next_batch(BatchDemand::new(2, 256 * 1024), &context)
            .expect("finished"),
        Backpressure::Finished
    ));
    assert_eq!(context.live_batches(), 0);
}

#[test]
fn demand_failure_restores_every_column_cursor() {
    let long_name = "A".repeat(80_000);
    let bytes = single_atom_file(&long_name);
    let context = ExecutionContext::default();
    let mut batches = BcifBatchSource::new(
        InputBuffer::from_bytes(bytes),
        ReadOptions::default(),
        DatasetId::new(1),
        ChunkId::new(0),
        LogicalRow::new(0),
        97,
        &context,
    )
    .expect("batch source");
    let error = batches
        .next_batch(BatchDemand::new(1, 66_000), &context)
        .expect_err("dictionary must exceed the first demand");
    assert!(matches!(error, StructureBatchError::DemandTooSmall { .. }));
    let ready = batches
        .next_batch(BatchDemand::new(1, 256 * 1024), &context)
        .expect("retry after rollback");
    let Backpressure::Ready(batch) = ready else {
        panic!("retry was not ready");
    };
    assert_eq!(batch.batch().positions(), &[[1.0, 2.0, 3.0]]);
    let atom = batch.batch().dictionary().resolve(batch.batch().atoms()[0]);
    assert_eq!(atom, Some(long_name.as_str()));
}

#[test]
fn oversized_demand_adapts_to_the_memory_still_available() {
    let bytes = repeated_atom_file(4_096);
    let probe = ExecutionContext::default();
    let probe_source = BcifBatchSource::new(
        InputBuffer::from_bytes(bytes.clone()),
        ReadOptions::default(),
        DatasetId::new(1),
        ChunkId::new(0),
        LogicalRow::new(0),
        97,
        &probe,
    )
    .expect("probe source");
    let source_bytes = probe.reserved_bytes();
    drop(probe_source);
    assert_eq!(probe.reserved_bytes(), 0);

    let budget = MemoryBudget::new(source_bytes + 200_000).expect("non-zero budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("bounded context");
    let mut source = BcifBatchSource::new(
        InputBuffer::from_bytes(bytes),
        ReadOptions::default(),
        DatasetId::new(1),
        ChunkId::new(0),
        LogicalRow::new(0),
        97,
        &context,
    )
    .expect("bounded source");
    let ready = source
        .next_batch(BatchDemand::new(4_096, 16 * 1024 * 1024), &context)
        .expect("adaptive pull");
    let Backpressure::Ready(batch) = ready else {
        panic!("a smaller batch should fit");
    };
    assert!(!batch.batch().is_empty());
    assert!(batch.batch().rows() < 4_096);
    assert!(context.reserved_bytes() <= budget.bytes());
    drop(batch);
    drop(source);
    assert_eq!(context.reserved_bytes(), 0);
}

fn atom_file() -> Vec<u8> {
    let rows = 3;
    let columns = vec![
        text("group_PDB", &["ATOM", "ATOM", "HETATM"]),
        integer("id", &[1, 2, 3]),
        text("type_symbol", &["C", "N", "O"]),
        text("label_atom_id", &["CA", "N", "O"]),
        nested_run_text("label_alt_id"),
        text("label_comp_id", &["GLY", "GLY", "HOH"]),
        text("label_asym_id", &["A", "A", "B"]),
        integer("label_seq_id", &[1, 1, 2]),
        float("Cartn_x", &[1.0, 4.0, 7.0]),
        float("Cartn_y", &[2.0, 5.0, 8.0]),
        float("Cartn_z", &[3.0, 6.0, 9.0]),
        float("occupancy", &[1.0, 1.0, 0.5]),
        float("B_iso_or_equiv", &[10.0, 11.0, 12.0]),
        integer("pdbx_PDB_model_num", &[1, 1, 1]),
    ];
    let file = EncodedFile {
        version: "0.3.0".to_owned(),
        encoder: "batch-test".to_owned(),
        data_blocks: vec![EncodedBlock {
            header: "test".to_owned(),
            categories: vec![EncodedCategory {
                name: "_atom_site".to_owned(),
                row_count: rows,
                columns,
            }],
        }],
    };
    crate::messagepack::consume_to_vec_named(file).expect("serialize fixture")
}

fn single_atom_file(atom_name: &str) -> Vec<u8> {
    let columns = vec![
        text("group_PDB", &["ATOM"]),
        integer("id", &[1]),
        text("type_symbol", &["C"]),
        text("label_atom_id", &[atom_name]),
        text("label_comp_id", &["GLY"]),
        text("label_asym_id", &["A"]),
        integer("label_seq_id", &[1]),
        float("Cartn_x", &[1.0]),
        float("Cartn_y", &[2.0]),
        float("Cartn_z", &[3.0]),
        integer("pdbx_PDB_model_num", &[1]),
    ];
    let file = EncodedFile {
        version: "0.3.0".to_owned(),
        encoder: "rollback-test".to_owned(),
        data_blocks: vec![EncodedBlock {
            header: "test".to_owned(),
            categories: vec![EncodedCategory {
                name: "_atom_site".to_owned(),
                row_count: 1,
                columns,
            }],
        }],
    };
    crate::messagepack::consume_to_vec_named(file).expect("serialize fixture")
}

fn repeated_atom_file(rows: usize) -> Vec<u8> {
    let atoms = vec!["CA"; rows];
    let groups = vec!["ATOM"; rows];
    let elements = vec!["C"; rows];
    let components = vec!["GLY"; rows];
    let chains = vec!["A"; rows];
    let integers = (0..rows)
        .map(|row| i64::try_from(row).expect("test row fits i64") + 1)
        .collect::<Vec<_>>();
    let coordinates = (0..rows)
        .map(|row| f64::from(u32::try_from(row).expect("test row fits u32")) + 1.0)
        .collect::<Vec<_>>();
    let columns = vec![
        text("group_PDB", &groups),
        integer("id", &integers),
        text("type_symbol", &elements),
        text("label_atom_id", &atoms),
        text("label_comp_id", &components),
        text("label_asym_id", &chains),
        integer("label_seq_id", &integers),
        float("Cartn_x", &coordinates),
        float("Cartn_y", &coordinates),
        float("Cartn_z", &coordinates),
        integer("pdbx_PDB_model_num", &vec![1; rows]),
    ];
    let file = EncodedFile {
        version: "0.3.0".to_owned(),
        encoder: "adaptive-test".to_owned(),
        data_blocks: vec![EncodedBlock {
            header: "test".to_owned(),
            categories: vec![EncodedCategory {
                name: "_atom_site".to_owned(),
                row_count: rows,
                columns,
            }],
        }],
    };
    crate::messagepack::consume_to_vec_named(file).expect("serialize fixture")
}

fn text(name: &str, values: &[&str]) -> EncodedColumn {
    let values = values
        .iter()
        .map(|value| (*value).to_owned())
        .collect::<Vec<_>>();
    column(name, encode_strings(&values).expect("string encoding"))
}

fn integer(name: &str, values: &[i64]) -> EncodedColumn {
    column(name, encode_integers(values).expect("integer encoding"))
}

fn float(name: &str, values: &[f64]) -> EncodedColumn {
    column(name, encode_floats(values).expect("float encoding"))
}

fn nested_run_text(name: &str) -> EncodedColumn {
    let offsets = encode_integers(&[0, 1]).expect("offset encoding");
    column(
        name,
        EncodedData {
            encoding: vec![Encoding::StringArray {
                data_encoding: vec![
                    Encoding::RunLength {
                        src_type: DataType::Int8,
                        src_size: 3,
                    },
                    Encoding::RunLength {
                        src_type: DataType::Int8,
                        src_size: 4,
                    },
                    Encoding::ByteArray {
                        r#type: DataType::Int8,
                    },
                ],
                string_data: "A".to_owned(),
                offset_encoding: offsets.encoding,
                offsets: offsets.data,
            }],
            data: [-1_i8, 1, 1, 1, 0, 1, 2, 1]
                .into_iter()
                .map(|value| value.to_le_bytes()[0])
                .collect(),
        },
    )
}

fn column(name: &str, data: EncodedData) -> EncodedColumn {
    EncodedColumn {
        name: name.to_owned(),
        data,
        mask: None,
    }
}
