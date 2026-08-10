//! Measures the safe snapshot constructor over a fixed local file.

use pdbiox_mmap::MappedFile;
use std::fs::OpenOptions;
use std::hint::black_box;
use std::io::{self, Write};
use std::time::Instant;

const SNAPSHOT_BYTES: usize = 8 * 1024 * 1024;
const ITERATIONS: usize = 10;

fn main() -> io::Result<()> {
    let path = std::env::temp_dir().join(format!("pdbiox-mmap-benchmark-{}", std::process::id()));
    let mut output = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(&path)?;
    output.write_all(&vec![17_u8; SNAPSHOT_BYTES])?;
    output.flush()?;

    let started = Instant::now();
    for _ in 0..ITERATIONS {
        black_box(MappedFile::new(&output)?);
    }
    let elapsed = started.elapsed();
    let total_bytes = SNAPSHOT_BYTES
        .checked_mul(ITERATIONS)
        .ok_or_else(|| io::Error::other("benchmark byte count overflow"))?;
    let mib = total_bytes / (1024 * 1024);
    eprintln!("snapshot: {mib} MiB in {elapsed:?}");
    std::fs::remove_file(path)?;
    Ok(())
}
