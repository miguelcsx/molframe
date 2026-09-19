use super::*;
use arrow::ffi_stream::ArrowArrayStreamReader;
use arrow::record_batch::RecordBatch;
use arrow::{datatypes::Schema, error::ArrowError};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Clone)]
struct ProbeSource {
    schema: SchemaRef,
    produced: Arc<AtomicUsize>,
    batches: usize,
    fail_at: Option<usize>,
}

impl ProbeSource {
    fn new(batches: usize, fail_at: Option<usize>) -> Self {
        Self {
            schema: Arc::new(Schema::empty()),
            produced: Arc::new(AtomicUsize::new(0)),
            batches,
            fail_at,
        }
    }
}

impl ArrowTableExport for ProbeSource {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn batch_count(&self) -> usize {
        self.batches
    }

    fn batch(&self, index: usize) -> Result<RecordBatch> {
        self.produced.fetch_add(1, Ordering::SeqCst);
        if self.fail_at == Some(index) {
            return Err(ArrowError::InvalidArgumentError(
                "deferred probe failure".to_owned(),
            ));
        }
        Ok(RecordBatch::new_empty(self.schema.clone()))
    }
}

#[test]
fn c_stream_materialises_exactly_one_batch_per_pull() {
    let source = ProbeSource::new(3, None);
    let produced = source.produced.clone();
    let stream = match source.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };
    assert_eq!(produced.load(Ordering::SeqCst), 0);

    let mut reader = match ArrowArrayStreamReader::try_new(stream.into_ffi()) {
        Ok(reader) => reader,
        Err(error) => panic!("stream import failed: {error}"),
    };
    assert_eq!(produced.load(Ordering::SeqCst), 0);

    for expected in 1..=3 {
        let Some(batch) = reader.next() else {
            panic!("batch {expected} absent")
        };
        assert!(batch.is_ok());
        assert_eq!(produced.load(Ordering::SeqCst), expected);
    }
    assert!(reader.next().is_none());
    assert_eq!(produced.load(Ordering::SeqCst), 3);
}

#[test]
fn c_stream_defers_chunk_errors_until_the_failing_pull() {
    let source = ProbeSource::new(3, Some(1));
    let produced = source.produced.clone();
    let stream = match source.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };
    let mut reader = match ArrowArrayStreamReader::try_new(stream.into_ffi()) {
        Ok(reader) => reader,
        Err(error) => panic!("stream import failed: {error}"),
    };
    assert_eq!(produced.load(Ordering::SeqCst), 0);
    assert!(reader.next().is_some_and(|batch| batch.is_ok()));
    assert_eq!(produced.load(Ordering::SeqCst), 1);
    let Some(error) = reader.next() else {
        panic!("deferred error absent")
    };
    assert!(error.is_err());
    assert_eq!(produced.load(Ordering::SeqCst), 2);
    assert!(reader.next().is_none());
    assert_eq!(produced.load(Ordering::SeqCst), 2);
}

#[test]
fn explicit_record_batch_export_materialises_every_batch_before_returning() {
    let source = ProbeSource::new(3, None);
    let produced = source.produced.clone();
    let batches = match source.record_batches() {
        Ok(batches) => batches,
        Err(error) => panic!("record batch export failed: {error}"),
    };
    assert_eq!(batches.len(), 3);
    assert_eq!(produced.load(Ordering::SeqCst), 3);
}

#[test]
fn c_stream_release_clears_callbacks_before_rust_drop() {
    let source = ProbeSource::new(1, None);
    let stream = match source.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };
    let mut ffi = stream.into_ffi();
    let Some(release) = ffi.release else {
        panic!("stream release callback absent")
    };
    assert!(ffi.get_schema.is_some());
    assert!(ffi.get_next.is_some());
    assert!(ffi.get_last_error.is_some());
    assert!(!ffi.private_data.is_null());

    // SAFETY: `ffi` is the live stream value produced above. This test invokes
    // its C release callback once and then verifies that Rust Drop is inert.
    unsafe { release(&raw mut ffi) };

    assert!(ffi.release.is_none());
    assert!(ffi.get_schema.is_none());
    assert!(ffi.get_next.is_none());
    assert!(ffi.get_last_error.is_none());
    drop(ffi);
}

#[test]
fn dropping_an_unconsumed_c_stream_does_not_materialise_batches() {
    let source = ProbeSource::new(4, None);
    let produced = Arc::clone(&source.produced);
    let stream = match source.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };

    drop(stream);

    assert_eq!(produced.load(Ordering::SeqCst), 0);
}
