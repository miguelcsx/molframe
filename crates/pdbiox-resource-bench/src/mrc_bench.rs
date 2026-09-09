//! Sparse-file resource probe for range-addressable density-map reads.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use super::{ResourceRecord, measure_case};

const HEADER_BYTES: usize = 1_024;
const AXIS: u32 = 1_024;
const AXIS_I32: i32 = 1_024;
const AXIS_F32: f32 = 1_024.0;
const BLOCK: usize = 96;

pub(super) fn run_mrc_block_1g() -> Result<ResourceRecord, String> {
    let fixture = SparseMrc::create()?;
    read_file("mrc_block_1g", &fixture.path)
}

pub(super) fn read_file(name: &'static str, path: &Path) -> Result<ResourceRecord, String> {
    measure_case(name, || {
        let mut reader =
            pdbiox::xtal::MrcBlockReader::open(path, pdbiox::xtal::MrcBlockOptions::default())
                .map_err(|error| format!("MRC block reader open failed: {error}"))?;
        let mut values = Vec::new();
        reader
            .read_block_into([113, 271, 509], [BLOCK; 3], &mut values)
            .map_err(|error| format!("first MRC block read failed: {error}"))?;
        let pointer = values.as_ptr();
        reader
            .read_block_into([701, 401, 17], [BLOCK; 3], &mut values)
            .map_err(|error| format!("second MRC block read failed: {error}"))?;
        if values.as_ptr() != pointer {
            return Err("MRC block output reallocated between equal requests".to_owned());
        }
        u64::try_from(values.len()).map_err(|_| "MRC block length exceeds u64".to_owned())
    })
}

struct SparseMrc {
    path: PathBuf,
}

impl SparseMrc {
    fn create() -> Result<Self, String> {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|error| format!("system clock cannot name MRC fixture: {error}"))?
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "pdbiox-mrc-resource-{}-{nonce}.mrc",
            std::process::id(),
        ));
        let mut file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&path)
            .map_err(|error| format!("cannot create sparse MRC fixture: {error}"))?;
        let mut header = [0_u8; HEADER_BYTES];
        write_i32_triplet(&mut header, 0, [AXIS_I32; 3]);
        write_i32(&mut header, 12, 0);
        write_i32_triplet(&mut header, 28, [AXIS_I32; 3]);
        write_f32_triplet(&mut header, 40, [AXIS_F32; 3]);
        write_f32_triplet(&mut header, 52, [90.0; 3]);
        write_i32_triplet(&mut header, 64, [1, 2, 3]);
        write_i32(&mut header, 88, 1);
        header[208..212].copy_from_slice(b"MAP ");
        header[212..216].copy_from_slice(&[0x44, 0x44, 0, 0]);
        file.write_all(&header)
            .map_err(|error| format!("cannot write sparse MRC header: {error}"))?;
        let values = u64::from(AXIS)
            .checked_pow(3)
            .ok_or_else(|| "sparse MRC size overflow".to_owned())?;
        let length = u64::try_from(HEADER_BYTES)
            .ok()
            .and_then(|header| header.checked_add(values))
            .ok_or_else(|| "sparse MRC length overflow".to_owned())?;
        file.set_len(length)
            .map_err(|error| format!("cannot size sparse MRC fixture: {error}"))?;
        drop(file);
        Ok(Self { path })
    }
}

impl Drop for SparseMrc {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

fn write_i32(bytes: &mut [u8], offset: usize, value: i32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn write_i32_triplet(bytes: &mut [u8], offset: usize, values: [i32; 3]) {
    for (index, value) in values.into_iter().enumerate() {
        write_i32(bytes, offset + index * 4, value);
    }
}

fn write_f32_triplet(bytes: &mut [u8], offset: usize, values: [f32; 3]) {
    for (index, value) in values.into_iter().enumerate() {
        bytes[offset + index * 4..offset + index * 4 + 4].copy_from_slice(&value.to_le_bytes());
    }
}
