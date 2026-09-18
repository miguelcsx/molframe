use std::io::{Cursor, Read, Seek, SeekFrom};

use molframe_core::structure::UnitCell;

use crate::DensityMap;

use super::*;

fn map(dimensions: [usize; 3]) -> DensityMap {
    let count = dimensions.into_iter().product();
    DensityMap {
        dimensions,
        starts: [0; 3],
        sampling: dimensions,
        cell: UnitCell {
            lengths: dimensions
                .map(|value| f64::from(u32::try_from(value).expect("test dimension fits in u32"))),
            angles: [90.0; 3],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: vec!["brick provider".into()],
        extended_header: Vec::new(),
        values: (0..count)
            .map(|value| {
                let bounded = value % (usize::from(u16::MAX) + 1);
                f32::from(u16::try_from(bounded).expect("bounded test value fits in u16"))
            })
            .collect(),
    }
}

fn provider(
    dimensions: [usize; 3],
    first_chunk: u64,
    first_brick: u64,
) -> MrcBrickProvider<Cursor<Vec<u8>>> {
    let bytes = map(dimensions).to_mrc_bytes().expect("test map encodes");
    let reader = MrcBlockReader::new(Cursor::new(bytes), MrcBlockOptions::default())
        .expect("test reader opens");
    MrcBrickProvider::new(
        reader,
        DatasetId::new(41),
        ChunkId::new(first_chunk),
        MapBrickId::new(first_brick),
        MrcBrickOptions {
            interior: [3, 3, 2],
            halo: 1,
            generation: 7,
            budget: MrcBrickBudget::default(),
        },
    )
    .expect("test provider builds")
}

#[test]
fn boundary_halos_repeat_the_nearest_source_voxel() {
    let mut provider = provider([5, 4, 3], 10, 20);
    let payload = provider
        .read_brick(MapBrickId::new(20))
        .expect("first brick reads");

    assert_eq!(payload.descriptor().metadata.shape.stored, [5, 5, 4]);
    assert_eq!(payload.values()[0].to_bits(), 0.0_f32.to_bits());
    assert_eq!(payload.values()[1].to_bits(), 0.0_f32.to_bits());
    assert_eq!(payload.values()[5].to_bits(), 0.0_f32.to_bits());
    assert_eq!(payload.values()[25].to_bits(), 0.0_f32.to_bits());
    assert_eq!(
        payload.actual_range().map(f32::to_bits),
        [0.0_f32.to_bits(), 58.0_f32.to_bits()]
    );
}

#[test]
fn positive_tail_bricks_keep_small_interiors_and_full_halos() {
    let provider = provider([5, 4, 3], 10, 20);
    assert_eq!(provider.brick_count(), 8);
    let descriptor = provider
        .descriptor(MapBrickId::new(27))
        .expect("tail descriptor exists");

    assert_eq!(descriptor.metadata.address.origin, [3, 3, 2]);
    assert_eq!(descriptor.metadata.shape.interior, [2, 1, 1]);
    assert_eq!(descriptor.metadata.shape.stored, [4, 3, 3]);
    assert_eq!(descriptor.metadata.shape.voxel_count, 36);
}

#[test]
fn global_brick_and_chunk_ids_are_not_narrowed_to_u32() {
    let first = u64::from(u32::MAX) + 90;
    let provider = provider([5, 4, 3], first, first + 100);
    let descriptor = provider
        .descriptor(MapBrickId::new(first + 107))
        .expect("large global identity resolves");

    assert_eq!(descriptor.chunk.get(), first + 7);
    assert_eq!(descriptor.metadata.id.get(), first + 107);
    assert_eq!(descriptor.metadata.address.origin, [3, 3, 2]);
}

#[derive(Debug)]
struct CountingSource {
    cursor: Cursor<Vec<u8>>,
    bytes_read: usize,
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
fn one_request_reads_and_retains_only_its_working_set() {
    let bytes = map([64, 64, 64])
        .to_mrc_bytes()
        .expect("bounded test map encodes");
    let complete_bytes = bytes.len();
    let source = CountingSource {
        cursor: Cursor::new(bytes),
        bytes_read: 0,
    };
    let reader =
        MrcBlockReader::new(source, MrcBlockOptions::default()).expect("counted reader opens");
    let mut provider = MrcBrickProvider::new(
        reader,
        DatasetId::new(1),
        ChunkId::new(2),
        MapBrickId::new(3),
        MrcBrickOptions {
            interior: [8; 3],
            halo: 1,
            generation: 0,
            budget: MrcBrickBudget {
                max_payload_bytes: 4_000,
                max_working_set_bytes: 5_000,
            },
        },
    )
    .expect("bounded provider builds");
    let payload = provider
        .read_brick(MapBrickId::new(3))
        .expect("bounded brick reads");
    let source = provider.into_reader().into_inner();

    assert_eq!(payload.payload_bytes(), 4_000);
    assert_eq!(payload.values().len(), 1_000);
    assert!(source.bytes_read < complete_bytes / 20);
}

#[test]
fn payload_budget_is_enforced_before_any_brick_request() {
    let bytes = map([8, 8, 8]).to_mrc_bytes().expect("test map encodes");
    let reader =
        MrcBlockReader::new(Cursor::new(bytes), MrcBlockOptions::default()).expect("reader opens");
    let result = MrcBrickProvider::new(
        reader,
        DatasetId::new(1),
        ChunkId::new(2),
        MapBrickId::new(3),
        MrcBrickOptions {
            interior: [8; 3],
            halo: 1,
            generation: 0,
            budget: MrcBrickBudget {
                max_payload_bytes: 3_999,
                max_working_set_bytes: 5_000,
            },
        },
    );
    assert!(matches!(result, Err(MrcBrickError::PayloadBudget { .. })));
}
