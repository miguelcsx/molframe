use super::*;
use std::fs::OpenOptions;
use std::io::{Read, Seek, Write};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

const LARGE_BYTES: usize = 8 * 1024 * 1024;

fn temporary_file(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("molframe-map-{name}-{}", std::process::id()))
}

#[test]
fn safe_new_is_an_owned_snapshot() {
    let path = temporary_file("safe-new");
    if let Err(error) = std::fs::write(&path, b"mapped bytes") {
        panic!("write failed: {error}")
    }
    let input = match File::open(&path) {
        Ok(file) => file,
        Err(error) => panic!("open failed: {error}"),
    };
    let snapshot = match MappedFile::new(&input) {
        Ok(snapshot) => snapshot,
        Err(error) => panic!("snapshot failed: {error}"),
    };
    assert!(matches!(&snapshot.mapping, Mapping::Snapshot(_)));

    drop(input);
    if let Err(error) = std::fs::write(&path, b"changed") {
        panic!("rewrite failed: {error}")
    }
    assert_eq!(snapshot.as_bytes(), b"mapped bytes");

    drop(snapshot);
    if let Err(error) = std::fs::remove_file(path) {
        panic!("cleanup failed: {error}")
    }
}

#[test]
fn large_snapshot_and_unchecked_mapping_have_distinct_backing() {
    let path = temporary_file("large");
    let pattern: Vec<u8> = (0_u8..=250).collect();
    let expected: Vec<u8> = pattern.iter().copied().cycle().take(LARGE_BYTES).collect();
    if let Err(error) = std::fs::write(&path, &expected) {
        panic!("write failed: {error}")
    }
    let input = match File::open(&path) {
        Ok(file) => file,
        Err(error) => panic!("open failed: {error}"),
    };
    let snapshot = match MappedFile::snapshot(&input) {
        Ok(snapshot) => snapshot,
        Err(error) => panic!("snapshot failed: {error}"),
    };
    // SAFETY: this test owns the path and does not modify the file until after
    // the file-backed mapping has been dropped.
    let mapped = match unsafe { MappedFile::map_file_unchecked(&input) } {
        Ok(mapped) => mapped,
        Err(error) => panic!("unchecked mapping failed: {error}"),
    };

    assert!(matches!(&snapshot.mapping, Mapping::Snapshot(_)));
    assert!(matches!(&mapped.mapping, Mapping::FileBacked { .. }));
    assert_eq!(snapshot.as_bytes(), expected);
    assert_eq!(mapped.as_bytes(), expected);

    drop(input);
    assert_eq!(&mapped[..pattern.len()], pattern);
    drop(mapped);

    let mut output = match OpenOptions::new().write(true).open(&path) {
        Ok(file) => file,
        Err(error) => panic!("reopen failed: {error}"),
    };
    if let Err(error) = output.rewind() {
        panic!("rewind failed: {error}")
    }
    if let Err(error) = output.write_all(b"changed") {
        panic!("rewrite failed: {error}")
    }
    assert_eq!(snapshot.as_bytes(), expected);

    drop(snapshot);
    drop(output);
    if let Err(error) = std::fs::remove_file(path) {
        panic!("cleanup failed: {error}")
    }
}

#[test]
fn safe_new_handles_an_empty_file_without_mapping() {
    let path = temporary_file("empty");
    if let Err(error) = std::fs::write(&path, []) {
        panic!("write failed: {error}")
    }
    let input = match File::open(&path) {
        Ok(file) => file,
        Err(error) => panic!("open failed: {error}"),
    };
    let mapped = match MappedFile::new(&input) {
        Ok(mapped) => mapped,
        Err(error) => panic!("snapshot failed: {error}"),
    };
    assert!(mapped.is_empty());
    assert!(matches!(&mapped.mapping, Mapping::Empty));
    drop(mapped);
    if let Err(error) = std::fs::remove_file(path) {
        panic!("cleanup failed: {error}")
    }
}

#[test]
fn private_reader_snapshot_retries_interrupts_and_bounds_each_read() {
    let expected = (0_u8..=250)
        .cycle()
        .take(3 * 64 * 1024 + 17)
        .collect::<Vec<_>>();
    let largest_request = Arc::new(AtomicUsize::new(0));
    let reader = InterruptingReader {
        bytes: expected.clone(),
        position: 0,
        interrupted: false,
        largest_request: Arc::clone(&largest_request),
    };

    let snapshot = match MappedFile::private_snapshot_from_reader(reader) {
        Ok(snapshot) => snapshot,
        Err(error) => panic!("private snapshot failed: {error}"),
    };

    assert_eq!(snapshot.as_bytes(), expected);
    assert!(matches!(&snapshot.mapping, Mapping::FileBacked { .. }));
    assert_eq!(largest_request.load(Ordering::SeqCst), 64 * 1024);
}

#[test]
fn anonymous_snapshot_does_not_change_the_callers_cursor() {
    let path = temporary_file("cursor");
    if let Err(error) = std::fs::write(&path, b"cursor remains here") {
        panic!("write failed: {error}")
    }
    let mut input = match OpenOptions::new().read(true).open(&path) {
        Ok(file) => file,
        Err(error) => panic!("open failed: {error}"),
    };
    if let Err(error) = input.seek(std::io::SeekFrom::Start(7)) {
        panic!("seek failed: {error}")
    }

    let snapshot = match MappedFile::snapshot(&input) {
        Ok(snapshot) => snapshot,
        Err(error) => panic!("snapshot failed: {error}"),
    };
    let cursor = match input.stream_position() {
        Ok(cursor) => cursor,
        Err(error) => panic!("cursor query failed: {error}"),
    };

    assert_eq!(cursor, 7);
    assert_eq!(snapshot.as_bytes(), b"cursor remains here");
    drop(snapshot);
    drop(input);
    if let Err(error) = std::fs::remove_file(path) {
        panic!("cleanup failed: {error}")
    }
}

#[test]
fn unchecked_empty_file_does_not_create_an_os_mapping() {
    let path = temporary_file("unchecked-empty");
    if let Err(error) = std::fs::write(&path, []) {
        panic!("write failed: {error}")
    }
    let input = match File::open(&path) {
        Ok(file) => file,
        Err(error) => panic!("open failed: {error}"),
    };
    // SAFETY: the empty fixture is owned by this test and is not modified
    // during the returned value's lifetime.
    let mapped = match unsafe { MappedFile::map_file_unchecked(&input) } {
        Ok(mapped) => mapped,
        Err(error) => panic!("unchecked mapping failed: {error}"),
    };

    assert!(mapped.is_empty());
    assert!(matches!(&mapped.mapping, Mapping::Empty));
    drop(mapped);
    drop(input);
    if let Err(error) = std::fs::remove_file(path) {
        panic!("cleanup failed: {error}")
    }
}

struct InterruptingReader {
    bytes: Vec<u8>,
    position: usize,
    interrupted: bool,
    largest_request: Arc<AtomicUsize>,
}

impl Read for InterruptingReader {
    fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
        self.largest_request
            .fetch_max(destination.len(), Ordering::SeqCst);
        if !self.interrupted {
            self.interrupted = true;
            return Err(io::Error::from(io::ErrorKind::Interrupted));
        }
        let remaining = &self.bytes[self.position..];
        let count = remaining.len().min(destination.len());
        destination[..count].copy_from_slice(&remaining[..count]);
        self.position += count;
        Ok(count)
    }
}
