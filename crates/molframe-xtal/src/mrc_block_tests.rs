use std::io::{Cursor, Read, Seek, SeekFrom};

use crate::DensityMap;

use super::*;

fn map(dimensions: [usize; 3]) -> DensityMap {
    let count = dimensions.into_iter().product();
    DensityMap {
        dimensions,
        starts: [-2, 3, 7],
        sampling: dimensions,
        cell: UnitCell {
            lengths: dimensions
                .map(|value| f64::from(u32::try_from(value).expect("test dimension fits in u32"))),
            angles: [90.0; 3],
        },
        origin: [1.0, 2.0, 3.0],
        space_group: 1,
        labels: vec!["block reader".into()],
        extended_header: b"extended".to_vec(),
        values: (0..count)
            .map(|index| f32::from(u16::try_from(index).expect("test value fits in u16")) * 0.25)
            .collect(),
    }
}

fn expected_block(map: &DensityMap, origin: [usize; 3], dimensions: [usize; 3]) -> Vec<f32> {
    let mut values = Vec::new();
    for z in 0..dimensions[2] {
        for y in 0..dimensions[1] {
            for x in 0..dimensions[0] {
                let Some(value) = map.value([origin[0] + x, origin[1] + y, origin[2] + z]) else {
                    panic!("test block must remain inside the map");
                };
                values.push(value);
            }
        }
    }
    values
}

#[derive(Debug)]
struct CountingSource {
    cursor: Cursor<Vec<u8>>,
    bytes_read: usize,
}

impl CountingSource {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            cursor: Cursor::new(bytes),
            bytes_read: 0,
        }
    }
}

impl Read for CountingSource {
    fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
        let count = self.cursor.read(output)?;
        self.bytes_read += count;
        Ok(count)
    }
}

impl Seek for CountingSource {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
        self.cursor.seek(position)
    }
}

#[test]
fn construction_reads_only_the_fixed_header() {
    let bytes = map([32, 32, 32]).to_mrc_bytes().expect("map should encode");
    let source = CountingSource::new(bytes);
    let reader =
        MrcBlockReader::new(source, MrcBlockOptions::default()).expect("header should decode");
    assert_eq!(reader.descriptor().dimensions, [32, 32, 32]);
    assert_eq!(reader.into_inner().bytes_read, HEADER_BYTES);
}

#[test]
fn a_small_block_does_not_read_the_complete_map() {
    let bytes = map([32, 32, 32]).to_mrc_bytes().expect("map should encode");
    let total_bytes = bytes.len();
    let source = CountingSource::new(bytes);
    let mut reader =
        MrcBlockReader::new(source, MrcBlockOptions::default()).expect("header should decode");
    let mut output = Vec::new();
    reader
        .read_block_into([3, 5, 7], [4, 4, 4], &mut output)
        .expect("block should decode");
    let source = reader.into_inner();
    assert_eq!(output.len(), 64);
    assert!(source.bytes_read < total_bytes / 10);
}

#[test]
fn canonical_blocks_match_the_full_decoder_and_reuse_output_capacity() {
    let expected_map = map([9, 8, 7]);
    let bytes = expected_map.to_mrc_bytes().expect("map should encode");
    let mut reader = MrcBlockReader::new(Cursor::new(bytes), MrcBlockOptions::default())
        .expect("header should decode");
    assert_eq!(reader.descriptor().starts, expected_map.starts);
    assert_eq!(reader.descriptor().labels, expected_map.labels);
    assert_eq!(
        reader.descriptor().extended_header_bytes,
        expected_map.extended_header.len()
    );
    let origin = [2, 3, 1];
    let dimensions = [4, 3, 5];
    let mut output = Vec::new();
    reader
        .read_block_into(origin, dimensions, &mut output)
        .expect("first block should decode");
    assert_eq!(output, expected_block(&expected_map, origin, dimensions));
    let pointer = output.as_ptr();
    let capacity = output.capacity();
    reader
        .read_block_into([1, 2, 0], dimensions, &mut output)
        .expect("second block should decode");
    assert_eq!(output.as_ptr(), pointer);
    assert_eq!(output.capacity(), capacity);
}

#[test]
fn permuted_file_axes_return_canonical_x_fastest_blocks() {
    let mut bytes = map([7, 6, 5]).to_mrc_bytes().expect("map should encode");
    bytes[64..68].copy_from_slice(&2_i32.to_le_bytes());
    bytes[68..72].copy_from_slice(&1_i32.to_le_bytes());
    bytes[72..76].copy_from_slice(&3_i32.to_le_bytes());
    let full = DensityMap::from_mrc_bytes(&bytes).expect("full map should decode");
    let mut reader = MrcBlockReader::new(Cursor::new(bytes), MrcBlockOptions::default())
        .expect("header should decode");
    assert_eq!(reader.descriptor().dimensions, full.dimensions);
    assert_eq!(reader.descriptor().starts, full.starts);
    let origin = [1, 2, 1];
    let dimensions = [3, 4, 2];
    let mut output = Vec::new();
    reader
        .read_block_into(origin, dimensions, &mut output)
        .expect("permuted block should decode");
    assert_eq!(output, expected_block(&full, origin, dimensions));
}

#[test]
fn packed_nibbles_work_across_odd_stored_rows() {
    let source_map = map([5, 3, 2]);
    let count = source_map.values.len();
    let mut bytes = source_map.to_mrc_bytes().expect("map should encode");
    bytes[12..16].copy_from_slice(&101_i32.to_le_bytes());
    let data_offset = HEADER_BYTES + source_map.extended_header.len();
    bytes.truncate(data_offset);
    for pair in 0..count.div_ceil(2) {
        let low = (pair * 2) % 16;
        let high = (pair * 2 + 1) % 16;
        bytes.push(u8::try_from(low | (high << 4)).expect("packed test value fits in u8"));
    }
    let full = DensityMap::from_mrc_bytes(&bytes).expect("packed map should decode");
    let mut reader = MrcBlockReader::new(Cursor::new(bytes), MrcBlockOptions::default())
        .expect("packed header should decode");
    let origin = [1, 0, 0];
    let dimensions = [3, 3, 2];
    let mut output = Vec::new();
    reader
        .read_block_into(origin, dimensions, &mut output)
        .expect("packed block should decode");
    assert_eq!(output, expected_block(&full, origin, dimensions));
}

#[test]
fn invalid_regions_and_memory_budgets_fail_before_payload_reads() {
    let bytes = map([8, 8, 8]).to_mrc_bytes().expect("map should encode");
    let source = CountingSource::new(bytes.clone());
    let options = MrcBlockOptions::default().with_memory_limit(255);
    let mut reader = MrcBlockReader::new(source, options).expect("header should decode");
    let mut output = vec![42.0];
    assert_eq!(
        reader.read_block_into([0, 0, 0], [4, 4, 4], &mut output),
        Err(MrcError::MemoryLimit {
            required: 272,
            limit: 255,
        })
    );
    assert_eq!(output, [42.0]);
    assert_eq!(reader.into_inner().bytes_read, HEADER_BYTES);

    let mut reader = MrcBlockReader::new(Cursor::new(bytes.clone()), MrcBlockOptions::default())
        .expect("header should decode");
    assert_eq!(
        reader.read_block_into([7, 0, 0], [2, 1, 1], &mut Vec::new()),
        Err(MrcError::InvalidRegion)
    );
    // The caller's budget is the bound, so a large budget is honoured rather
    // than refused against a fixed ceiling; only a zero budget is invalid, and
    // it is rejected before the header is parsed.
    let error = MrcBlockReader::new(
        Cursor::new(vec![0; HEADER_BYTES]),
        MrcBlockOptions::default().with_memory_limit(0),
    )
    .expect_err("a zero memory limit must be rejected before header parsing");
    assert_eq!(error, MrcError::InvalidMemoryLimit { requested: 0 });

    let accepted = MrcBlockReader::new(
        Cursor::new(bytes.clone()),
        MrcBlockOptions::default().with_memory_limit(8_000_000_000),
    );
    assert!(
        accepted.is_ok(),
        "a budget past the former 500 MB ceiling must be accepted"
    );
}

#[test]
fn truncated_requested_rows_are_reported_without_panicking() {
    let mut bytes = map([8, 8, 8]).to_mrc_bytes().expect("map should encode");
    bytes.truncate(bytes.len() - 8);
    let mut reader = MrcBlockReader::new(Cursor::new(bytes), MrcBlockOptions::default())
        .expect("header remains available");
    assert_eq!(
        reader.read_block_into([0, 0, 7], [8, 8, 1], &mut Vec::new()),
        Err(MrcError::Truncated)
    );
}
