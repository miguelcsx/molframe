use super::*;
use crate::{Frame, StreamingReader};

#[test]
fn streaming_rmsf_matches_resident_reference_and_retains_its_charge() {
    let first = [[0.0; 3]; 10];
    let second = [[2.0; 3]; 10];
    let expected = molframe_geom::rmsf(&[&first, &second]).expect("reference");
    let mut reader = StreamingReader::new(
        [
            Frame {
                positions: first.to_vec(),
            },
            Frame {
                positions: second.to_vec(),
            },
        ]
        .into_iter(),
        10,
    );
    let context = ExecutionContext::default();
    let result = rmsf_stream(&mut reader, &context, 4096).expect("rmsf");
    assert_eq!(result.as_slice(), expected.as_slice());
    assert_eq!(context.reserved_bytes(), 80);
    drop(result);
    assert_eq!(context.reserved_bytes(), 0);
}
