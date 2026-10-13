use super::*;

fn sample() -> FrameIndex {
    let stamp = FileStamp {
        len: 1_000,
        modified_ns: 77,
    };
    FrameIndex::new(FrameIndex::XTC, 12, stamp, vec![0, 100, 350, 999]).expect("index")
}

#[test]
fn bytes_round_trip() {
    let index = sample();
    let again = FrameIndex::from_bytes(&index.to_bytes()).expect("parse");
    assert_eq!(again, index);
    assert_eq!(again.len(), 4);
    assert!(again.describes(FrameIndex::XTC, 12, index.stamp()));
    assert!(!again.describes(FrameIndex::TRR, 12, index.stamp()));
    assert!(!again.describes(FrameIndex::XTC, 13, index.stamp()));
}

#[test]
fn file_round_trip_and_missing_file() {
    let dir = tempfile::tempdir().expect("dir");
    let path = dir.path().join("a.idx");
    let index = sample();
    index.write(&path).expect("write");
    assert_eq!(FrameIndex::read(&path).expect("read"), index);
    assert!(matches!(
        FrameIndex::read(&dir.path().join("none")),
        Err(FrameIndexError::Io(_))
    ));
}

#[test]
fn offsets_must_be_ordered_and_inside_the_file() {
    let stamp = FileStamp {
        len: 100,
        modified_ns: 0,
    };
    for bad in [vec![5], vec![0, 0], vec![0, 50, 40], vec![0, 100]] {
        assert_eq!(
            FrameIndex::new(FrameIndex::TRR, 3, stamp, bad),
            Err(FrameIndexError::Malformed)
        );
    }
}

#[test]
fn forged_count_is_refused_without_allocating() {
    let mut bytes = sample().to_bytes();
    bytes[40..48].copy_from_slice(&u64::MAX.to_le_bytes());
    let body = bytes.len() - 8;
    let sum = checksum(&bytes[..body]);
    bytes[body..].copy_from_slice(&sum.to_le_bytes());
    assert_eq!(
        FrameIndex::from_bytes(&bytes),
        Err(FrameIndexError::Malformed)
    );
}

#[test]
fn damaged_bytes_are_refused() {
    let good = sample().to_bytes();
    for cut in 0..good.len() {
        assert!(FrameIndex::from_bytes(&good[..cut]).is_err(), "cut {cut}");
    }
    for at in 0..good.len() {
        let mut bad = good.clone();
        bad[at] ^= 0x40;
        assert!(FrameIndex::from_bytes(&bad).is_err(), "flip {at}");
    }
}

#[test]
fn stamp_reads_a_real_file() {
    let dir = tempfile::tempdir().expect("dir");
    let path = dir.path().join("f");
    std::fs::write(&path, [1, 2, 3]).expect("write");
    let stamp = FileStamp::of(&path).expect("stamp");
    assert_eq!(stamp.len, 3);
    assert!(FileStamp::of(&dir.path().join("none")).is_err());
}
