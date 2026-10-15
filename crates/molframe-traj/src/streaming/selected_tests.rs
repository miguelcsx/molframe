use super::*;
use crate::{Frame, MemoryReader, StreamingReader, Trajectory};

const ATOMS: usize = 4;

fn small(value: usize) -> f32 {
    f32::from(u16::try_from(value).expect("small"))
}

fn frames(count: usize) -> Vec<Frame> {
    (0..count)
        .map(|frame| Frame {
            positions: (0..ATOMS)
                .map(|atom| [small(frame), small(atom), small(frame * 10 + atom)])
                .collect(),
        })
        .collect()
}

fn collect<R: TrajectoryReader>(reader: &mut R) -> Vec<(usize, Vec<[f32; 3]>)> {
    let mut timestep = Timestep::default();
    let mut out = Vec::new();
    while reader.read_next(&mut timestep).expect("read") {
        out.push((timestep.frame, timestep.positions.clone()));
    }
    out
}

#[test]
fn seekable_source_yields_the_strided_frames_with_source_ordinals() {
    let trajectory = Trajectory::from_frames(frames(10)).expect("trajectory");
    let slice = FrameSlice {
        start: 1,
        stop: 9,
        step: 3,
    };
    let mut reader =
        SelectedReader::new(MemoryReader::new(&trajectory), slice, None).expect("reader");
    assert_eq!(reader.n_frames(), Some(3));
    let got = collect(&mut reader);
    assert_eq!(got.iter().map(|f| f.0).collect::<Vec<_>>(), [1, 4, 7]);
    assert!(
        got[1].1[2]
            .iter()
            .zip([4, 2, 42])
            .all(|(a, b)| (a - small(b)).abs() < f32::EPSILON)
    );
}

#[test]
fn forward_only_source_gives_the_same_frames_by_decoding_and_discarding() {
    let slice = FrameSlice {
        start: 2,
        stop: usize::MAX,
        step: 4,
    };
    let mut reader = SelectedReader::new(
        StreamingReader::new(frames(11).into_iter(), ATOMS),
        slice,
        Some(vec![3, 0]),
    )
    .expect("reader");
    assert_eq!(reader.n_atoms(), 2);
    let got = collect(&mut reader);
    assert_eq!(got.iter().map(|f| f.0).collect::<Vec<_>>(), [2, 6, 10]);
    assert_eq!(got[0].1, vec![[2.0, 3.0, 23.0], [2.0, 0.0, 20.0]]);
}

#[test]
fn atom_subset_follows_the_requested_order_and_rejects_out_of_range_indices() {
    let trajectory = Trajectory::from_frames(frames(2)).expect("trajectory");
    let bad = SelectedReader::new(
        MemoryReader::new(&trajectory),
        FrameSlice::ALL,
        Some(vec![0, 4]),
    );
    assert!(matches!(
        bad.err(),
        Some(TrajectoryError::SelectionOutOfRange { index: 4, atoms: 4 })
    ));
    let zero = SelectedReader::new(
        MemoryReader::new(&trajectory),
        FrameSlice {
            step: 0,
            ..FrameSlice::ALL
        },
        None,
    );
    assert!(zero.is_err());
}

#[test]
fn dense_chunks_are_filled_in_place_and_split_at_the_requested_size() {
    let trajectory = Trajectory::from_frames(frames(7)).expect("trajectory");
    let mut reader = SelectedReader::new(
        MemoryReader::new(&trajectory),
        FrameSlice::ALL,
        Some(vec![1, 2]),
    )
    .expect("reader");
    let mut scratch = Timestep::default();
    let first = read_dense_chunk(&mut reader, &mut scratch, 3).expect("chunk");
    assert_eq!((first.frames(), first.atoms), (3, 2));
    assert_eq!(first.coordinates.len(), 3 * 2 * 3);
    assert_eq!(first.coordinates.capacity(), first.coordinates.len());
    assert!(first.times.iter().all(|time| time.is_nan()));
    assert!(first.cells.is_none());
    let second = read_dense_chunk(&mut reader, &mut scratch, 3).expect("chunk");
    let third = read_dense_chunk(&mut reader, &mut scratch, 3).expect("chunk");
    assert_eq!((second.indices, third.indices), (vec![3, 4, 5], vec![6]));
    assert_eq!(count_frames(&mut reader, &mut scratch), Ok(0));
}

#[test]
fn bounded_reads_refuse_a_subset_workspace_over_the_ceiling() {
    let mut reader = SelectedReader::new(
        StreamingReader::new(frames(2).into_iter(), ATOMS),
        FrameSlice::ALL,
        Some(vec![0]),
    )
    .expect("reader");
    let mut timestep = Timestep::default();
    assert!(matches!(
        reader.read_next_bounded(&mut timestep, 16),
        Err(TrajectoryError::MemoryLimit { .. })
    ));
}

#[test]
fn rewinding_a_selection_replays_it_for_seekable_sources_only() {
    let trajectory = Trajectory::from_frames(frames(6)).expect("trajectory");
    let slice = FrameSlice {
        start: 0,
        stop: 6,
        step: 2,
    };
    let mut reader =
        SelectedReader::new(MemoryReader::new(&trajectory), slice, None).expect("reader");
    assert_eq!(collect(&mut reader).len(), 3);
    reader.rewind().expect("rewind");
    assert_eq!(collect(&mut reader).len(), 3);
    let mut stream = SelectedReader::new(
        StreamingReader::new(frames(6).into_iter(), ATOMS),
        slice,
        None,
    )
    .expect("reader");
    assert_eq!(collect(&mut stream).len(), 3);
    assert_eq!(
        stream.rewind(),
        Err(TrajectoryError::RandomAccessUnavailable)
    );
}
