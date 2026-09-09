use super::{
    ChainedReader, MemoryReader, RandomAccess, StreamingReader, Timestep, TrajectoryError,
    TrajectoryReader,
};
use crate::{Frame, Trajectory};

fn frame(value: f32) -> Frame {
    Frame {
        positions: vec![[value, 0.0, 0.0]],
    }
}

#[test]
fn stream_is_forward_only_and_reuses_the_caller_buffer() {
    let frames = (0_u16..10_000).map(|index| frame(f32::from(index)));
    let mut reader = StreamingReader::new(frames, 1);
    let mut timestep = Timestep::default();
    let mut pointer = None;
    let mut count = 0usize;
    while reader
        .read_next(&mut timestep)
        .unwrap_or_else(|error| panic!("stream failed: {error}"))
    {
        let current = timestep.positions.as_ptr();
        if let Some(first) = pointer {
            assert_eq!(current, first, "coordinate allocation should be reused");
        }
        pointer = Some(current);
        count += 1;
    }
    assert_eq!(count, 10_000);
    assert_eq!(reader.n_frames(), None);
    assert_eq!(reader.random_access(), RandomAccess::None);
    assert_eq!(
        reader.seek(0),
        Err(TrajectoryError::RandomAccessUnavailable)
    );
}

#[test]
fn memory_source_declares_and_performs_random_access() {
    let trajectory =
        Trajectory::from_frames(vec![frame(0.0), frame(1.0)]).expect("fixed-width trajectory");
    let mut reader = MemoryReader::new(&trajectory);
    let mut timestep = Timestep::default();
    reader
        .seek(1)
        .unwrap_or_else(|error| panic!("seek failed: {error}"));
    assert!(reader.read_next(&mut timestep).is_ok_and(|read| read));
    assert_eq!(timestep.frame, 1);
    assert_eq!(reader.random_access(), RandomAccess::Full);
}

#[test]
fn chained_sources_are_validated_and_numbered_globally() {
    let first = Trajectory::from_frames(vec![frame(0.0)]).expect("fixed-width trajectory");
    let second =
        Trajectory::from_frames(vec![frame(1.0), frame(2.0)]).expect("fixed-width trajectory");
    let mut chain = ChainedReader::new(vec![MemoryReader::new(&first), MemoryReader::new(&second)])
        .unwrap_or_else(|error| panic!("chain failed: {error}"));
    let mut timestep = Timestep::default();
    for expected in 0..3 {
        assert!(chain.read_next(&mut timestep).is_ok_and(|read| read));
        assert_eq!(timestep.frame, expected);
    }
    assert!(chain.read_next(&mut timestep).is_ok_and(|read| !read));
}

#[test]
fn chain_refuses_misaligned_atom_counts() {
    let first = Trajectory::from_frames(vec![frame(0.0)]).expect("fixed-width trajectory");
    let second = Trajectory::from_frames(vec![Frame {
        positions: vec![[0.0; 3]; 2],
    }])
    .expect("fixed-width trajectory");
    let result = ChainedReader::new(vec![MemoryReader::new(&first), MemoryReader::new(&second)]);
    assert!(matches!(
        result,
        Err(TrajectoryError::AtomCountMismatch {
            expected: 1,
            found: 2
        })
    ));
}

#[test]
fn an_unrepresentable_chain_length_is_reported_as_unknown_instead_of_wrapping() {
    struct DeclaredReader(usize);
    impl TrajectoryReader for DeclaredReader {
        fn format(&self) -> &'static str {
            "declared"
        }
        fn n_atoms(&self) -> usize {
            0
        }
        fn n_frames(&self) -> Option<usize> {
            Some(self.0)
        }
        fn units(&self) -> super::Units {
            super::Units::CANONICAL
        }
        fn random_access(&self) -> RandomAccess {
            RandomAccess::None
        }
        fn read_next(&mut self, _: &mut Timestep) -> Result<bool, TrajectoryError> {
            Ok(false)
        }
        fn seek(&mut self, _: usize) -> Result<(), TrajectoryError> {
            Err(TrajectoryError::RandomAccessUnavailable)
        }
    }
    let chain = ChainedReader::new(vec![DeclaredReader(usize::MAX), DeclaredReader(1)])
        .expect("matching topology");
    assert_eq!(chain.n_frames(), None);
}
