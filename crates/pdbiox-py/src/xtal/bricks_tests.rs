use std::io::Cursor;

use numpy::{PyArrayMethods, PyUntypedArrayMethods};
use pdbiox::UnitCell;

use super::*;

fn payload() -> pdbiox::xtal::ScalarBrickPayload {
    let map = pdbiox::xtal::DensityMap {
        dimensions: [4, 4, 4],
        starts: [0; 3],
        sampling: [4; 3],
        cell: UnitCell {
            lengths: [4.0; 3],
            angles: [90.0; 3],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values: (0_u16..64).map(f32::from).collect(),
    };
    let bytes = map.to_mrc_bytes().expect("test map encodes");
    let reader = pdbiox::xtal::MrcBlockReader::new(
        Cursor::new(bytes),
        pdbiox::xtal::MrcBlockOptions::default(),
    )
    .expect("test reader opens");
    let mut provider = pdbiox::xtal::MrcBrickProvider::new(
        reader,
        pdbiox::DatasetId::new(1),
        pdbiox::ChunkId::new(2),
        pdbiox::xtal::MapBrickId::new(3),
        pdbiox::xtal::MrcBrickOptions {
            interior: [2; 3],
            halo: 1,
            generation: 4,
            budget: pdbiox::xtal::MrcBrickBudget::default(),
        },
    )
    .expect("test provider builds");
    provider
        .read_brick(pdbiox::xtal::MapBrickId::new(3))
        .expect("test payload reads")
}

#[test]
fn scalar_values_are_zero_copy_contiguous_read_only_and_owner_retained() {
    Python::initialize();
    Python::attach(|py| {
        if py.import("numpy").is_err() {
            return;
        }
        let payload = payload();
        let pointer = payload.values().as_ptr();
        let owner = Bound::new(py, PyScalarBrickPayload(payload))
            .expect("payload should become a Python owner");
        let values = PyScalarBrickPayload::values(&owner);

        assert_eq!(values.data().cast_const(), pointer);
        assert!(values.is_c_contiguous());
        assert!(values.try_readwrite().is_err());

        drop(owner);
        let readonly = values.readonly();
        let slice = readonly
            .as_slice()
            .expect("owner-retained values remain contiguous");
        assert_eq!(slice.len(), 64);
        assert_eq!(slice[0].to_bits(), 0.0_f32.to_bits());
    });
}
