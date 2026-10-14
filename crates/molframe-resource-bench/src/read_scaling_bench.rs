//! Allocation guard for reading a large synthetic structure.
//!
//! The parallel readers cut their input into pieces, so their allocation count
//! must grow with the number of pieces and rows' worth of buffers, never with a
//! per-value cost. A ceiling on a 200k-atom read catches a reader that starts
//! allocating per value or per atom again.

use std::hint::black_box;

use super::{ResourceRecord, measure_case};

const ATOMS: u64 = 200_000;

pub(super) fn mmcif() -> Result<ResourceRecord, String> {
    let bytes = molframe_bench::cif_bytes_with_atoms(ATOMS);
    measure_case("mmcif_read_200k", || read(&bytes, "synthetic.cif"))
}

pub(super) fn bcif() -> Result<ResourceRecord, String> {
    let text = molframe_bench::cif_bytes_with_atoms(ATOMS);
    let (structure, _) = molframe::read_bytes(
        text.to_vec(),
        Some("synthetic.cif"),
        &molframe::ReadOptions::new(),
    )
    .map_err(|findings| format!("synthetic mmCIF failed to read: {findings:?}"))?;
    let bytes = molframe::write_bcif(&structure)
        .map_err(|findings| format!("synthetic BinaryCIF failed to write: {findings:?}"))?;
    measure_case("bcif_read_200k", || read(&bytes, "synthetic.bcif"))
}

fn read(bytes: &[u8], name: &'static str) -> Result<u64, String> {
    let (structure, _) =
        molframe::read_bytes(bytes.to_vec(), Some(name), &molframe::ReadOptions::new())
            .map_err(|findings| format!("{name} read failed: {findings:?}"))?;
    Ok(black_box(u64::from(structure.atom_count())))
}
