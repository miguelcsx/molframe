//! Criterion coverage for the safe read-only mmap snapshot constructor.

use criterion::{Criterion, black_box};
use pdbiox_mmap::MappedFile;
use std::fs::OpenOptions;
use std::io::{self, Write};

const SNAPSHOT_BYTES: usize = 8 * 1024 * 1024;

fn bench_snapshot(c: &mut Criterion) -> io::Result<()> {
    let path = std::env::temp_dir().join(format!("pdbiox-mmap-criterion-{}", std::process::id()));
    let mut output = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(&path)?;
    output.write_all(&vec![17_u8; SNAPSHOT_BYTES])?;
    output.flush()?;
    c.bench_function("mmap_snapshot", |b| {
        b.iter(|| black_box(MappedFile::new(&output)));
    });
    std::fs::remove_file(path)
}

fn main() -> io::Result<()> {
    let mut criterion = Criterion::default().configure_from_args();
    bench_snapshot(&mut criterion)?;
    criterion.final_summary();
    Ok(())
}
