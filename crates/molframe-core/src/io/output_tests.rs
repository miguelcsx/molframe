use super::*;
use std::cell::Cell;
use std::io::{self, Write};
use std::rc::Rc;

const BYTES: &[u8] = b"data_test\n_entry.id test\n#\n";

struct CountingWriter {
    bytes: Rc<Cell<usize>>,
    writes: Rc<Cell<usize>>,
}

impl Write for CountingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes.set(self.bytes.get() + bytes.len());
        self.writes.set(self.writes.get() + 1);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn output_path(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("molframe-output-{}-{name}", std::process::id()))
}

fn write_and_read(path: &Path) -> Vec<u8> {
    if let Err(error) = write_output(path, BYTES) {
        panic!("write failed: {error}")
    }
    let actual = match std::fs::read(path) {
        Ok(actual) => actual,
        Err(error) => panic!("read failed: {error}"),
    };
    if let Err(error) = std::fs::remove_file(path) {
        panic!("cleanup failed: {error}")
    }
    actual
}

#[test]
fn plain_output_writes_exactly_the_supplied_bytes() {
    let actual = write_and_read(&output_path("plain.cif"));
    assert_eq!(actual, BYTES);
}

#[test]
fn the_default_is_explicit_and_a_zero_limit_is_refused() {
    assert_eq!(
        OutputOptions::default().memory_limit_bytes,
        DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES
    );
    let path = output_path("invalid.cif");
    let options = OutputOptions::default().with_memory_limit(0);
    assert!(OutputSink::create(&path, options).is_err());
    assert!(
        !path.exists(),
        "invalid policy must refuse before file creation"
    );
}

#[test]
fn a_limit_larger_than_the_default_is_accepted() {
    let path = output_path("large-limit.cif");
    let options =
        OutputOptions::default().with_memory_limit(DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES * 100);
    assert!(OutputSink::create(&path, options).is_ok());
}

#[test]
fn bounded_sink_forwards_chunks_before_finish() {
    let written_bytes = Rc::new(Cell::new(0));
    let write_calls = Rc::new(Cell::new(0));
    let writer = CountingWriter {
        bytes: Rc::clone(&written_bytes),
        writes: Rc::clone(&write_calls),
    };
    let options = OutputOptions::default().with_memory_limit(1_024);
    let mut output =
        match OutputSink::from_writer(writer, OutputCompression::None, options, "counter".into()) {
            Ok(output) => output,
            Err(error) => panic!("sink construction failed: {error}"),
        };
    assert_eq!(output.buffer_capacity(), Some(1_024));
    for _ in 0..3 {
        if let Err(error) = output.write_all(&[7; 1_024]) {
            panic!("stream write failed: {error}");
        }
    }
    assert!(write_calls.get() > 0, "buffer must drain before finish");
    assert!(
        written_bytes.get() > 0,
        "bytes must reach the destination early"
    );
    if let Err(error) = output.finish() {
        panic!("stream finish failed: {error}");
    }
    assert_eq!(written_bytes.get(), 3_072);
}

#[cfg(feature = "gzip")]
#[test]
fn gzip_output_is_byte_deterministic_and_round_trips() {
    use std::io::Read;

    let first = write_and_read(&output_path("first.cif.gz"));
    let second = write_and_read(&output_path("second.cif.gz"));
    assert_eq!(first, second);

    let mut round_tripped = Vec::new();
    let mut decoder = flate2::read::GzDecoder::new(first.as_slice());
    if let Err(error) = decoder.read_to_end(&mut round_tripped) {
        panic!("gzip decode failed: {error}")
    }
    assert_eq!(round_tripped, BYTES);
}

#[cfg(feature = "zstd")]
#[test]
fn zstd_output_is_byte_deterministic_and_round_trips() {
    let first = write_and_read(&output_path("first.cif.zst"));
    let second = write_and_read(&output_path("second.cif.zst"));
    assert_eq!(first, second);

    let decoded = match zstd::stream::decode_all(first.as_slice()) {
        Ok(decoded) => decoded,
        Err(error) => panic!("zstandard decode failed: {error}"),
    };
    assert_eq!(decoded, BYTES);
}
