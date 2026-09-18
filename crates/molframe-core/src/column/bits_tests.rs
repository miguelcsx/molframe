use super::*;

#[test]
fn a_bit_set_reports_the_values_that_were_pushed_into_it() {
    let bits = BitVec::try_from_iter([true, false, true, true]).expect("small bit vector");
    assert_eq!(bits.len(), 4);
    assert_eq!(bits.get(0), Some(true));
    assert_eq!(bits.get(1), Some(false));
    assert_eq!(bits.get(4), None);
    assert_eq!(bits.count_ones(), 3);
}

#[test]
fn iterating_yields_set_positions_in_ascending_order() {
    let mut bits = BitVec::repeat(false, 200);
    for position in [0, 63, 64, 65, 199] {
        bits.set(position, true);
    }
    assert_eq!(bits.ones().collect::<Vec<_>>(), [0, 63, 64, 65, 199]);
}

#[test]
fn bits_past_the_end_never_show_up_in_a_count_after_inversion() {
    let mut bits = BitVec::repeat(false, 3);
    bits.invert();
    assert_eq!(bits.count_ones(), 3);
    assert!(bits.all());
    assert_eq!(bits.ones().collect::<Vec<_>>(), [0, 1, 2]);
}

#[test]
fn intersection_keeps_only_shared_bits_and_union_adds_them() {
    let mut left = BitVec::repeat(false, 100);
    let mut right = BitVec::repeat(false, 100);
    left.set(1, true);
    left.set(70, true);
    right.set(70, true);
    right.set(99, true);

    let mut both = left.clone();
    both.intersect_with(&right);
    assert_eq!(both.ones().collect::<Vec<_>>(), [70]);

    left.union_with(&right);
    assert_eq!(left.ones().collect::<Vec<_>>(), [1, 70, 99]);
}

#[test]
fn a_set_of_all_zeroes_reports_none_and_never_all() {
    let bits = BitVec::repeat(false, 10);
    assert!(bits.none());
    assert!(!bits.all());
    assert_eq!(bits.ones().next(), None);
}

#[test]
fn width_is_the_smallest_that_represents_the_largest_value() {
    assert_eq!(bit_width(0), 1);
    assert_eq!(bit_width(1), 1);
    assert_eq!(bit_width(2), 2);
    assert_eq!(bit_width(127), 7);
    assert_eq!(bit_width(128), 8);
    assert_eq!(bit_width(u64::MAX), 64);
}

#[test]
fn packed_values_survive_a_round_trip_at_every_width_they_fit() {
    let values: Vec<u64> = (0u64..37).map(|i| i * 3 % 100).collect();
    for width in [7, 8, 13, 32, 64] {
        let packed = pack(&values, width).expect("small packed buffer");
        for (index, expected) in values.iter().enumerate() {
            let read = unpack_one(&packed, width, u32::try_from(index).expect("small index"));
            assert_eq!(read, Some(*expected), "width {width} index {index}");
        }
    }
}

#[test]
fn reading_past_the_end_of_a_packed_buffer_yields_nothing() {
    let packed = pack(&[1, 2, 3], 8).expect("small packed buffer");
    assert_eq!(unpack_one(&packed, 8, 2), Some(3));
    assert_eq!(unpack_one(&packed, 8, 3), None);
}

#[test]
fn every_packed_value_reads_back_unchanged_at_every_width() {
    // Every width from 1 to 64 exercises a different byte span and bit offset,
    // including the nine-byte straddle that only a 64-bit value can reach.
    for width in 1..=64u8 {
        let mask = if width == 64 {
            u64::MAX
        } else {
            (1u64 << width) - 1
        };
        let values: Vec<u64> = (0..37u64)
            .map(|index| index.wrapping_mul(0x9E37_79B9_7F4A_7C15) & mask)
            .collect();

        let packed = pack(&values, width).expect("values pack");

        for (index, expected) in values.iter().enumerate() {
            let index = u32::try_from(index).expect("test index fits in u32");
            assert_eq!(
                unpack_one(&packed, width, index),
                Some(*expected),
                "width {width}, index {index}"
            );
        }
    }
}
