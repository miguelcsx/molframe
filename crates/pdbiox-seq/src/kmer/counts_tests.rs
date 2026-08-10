use super::{kmer_counts, minimizers};

#[test]
fn overlapping_windows_are_all_counted() {
    let counts = kmer_counts(b"AAAA", 2);
    assert_eq!(counts, vec![(b"AA".to_vec(), 3)]);
}

#[test]
fn distinct_words_come_back_sorted() {
    let counts = kmer_counts(b"ACGT", 1);
    let words: Vec<Vec<u8>> = counts.iter().map(|(word, _)| word.clone()).collect();
    assert_eq!(
        words,
        vec![b"A".to_vec(), b"C".to_vec(), b"G".to_vec(), b"T".to_vec()]
    );
    assert!(counts.iter().all(|(_, count)| *count == 1));
}

#[test]
fn a_repeated_word_accumulates() {
    // ACAC has windows AC, CA, AC — so AC appears twice.
    let counts = kmer_counts(b"ACAC", 2);
    let ac = counts.iter().find(|(word, _)| word == b"AC");
    let Some((_, count)) = ac else {
        panic!("expected an AC k-mer");
    };
    assert_eq!(*count, 2);
}

#[test]
fn a_window_longer_than_the_sequence_has_no_kmers() {
    assert!(kmer_counts(b"AC", 3).is_empty());
}

#[test]
fn a_zero_length_window_has_no_kmers() {
    assert!(kmer_counts(b"ACGT", 0).is_empty());
}

#[test]
fn minimizers_pick_the_smallest_kmer_of_each_window() {
    // k-mers of ACACAC: AC CA AC CA AC. In every 2-mer window the smallest is
    // "AC", so the minimizers are the AC positions 0, 2 and 4.
    let picked = minimizers(b"ACACAC", 2, 2);
    assert_eq!(
        picked,
        vec![
            (0, b"AC".to_vec()),
            (2, b"AC".to_vec()),
            (4, b"AC".to_vec()),
        ]
    );
}

#[test]
fn minimizers_break_ties_toward_the_left() {
    // k-mers of ACAC: AC CA AC. One window of all three ties two "AC"s; the
    // earlier position wins.
    let picked = minimizers(b"ACAC", 2, 3);
    assert_eq!(picked, vec![(0, b"AC".to_vec())]);
}

#[test]
fn a_sequence_too_short_for_a_window_has_no_minimizers() {
    assert!(minimizers(b"AC", 2, 2).is_empty());
}
