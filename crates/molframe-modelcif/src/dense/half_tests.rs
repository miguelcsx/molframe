use super::{f16_bits_to_f32, f32_to_f16_bits};

#[test]
fn every_finite_binary16_pattern_round_trips_through_single_precision() {
    for bits in 0..=u16::MAX {
        let value = f16_bits_to_f32(bits);
        if value.is_nan() {
            assert!(f16_bits_to_f32(f32_to_f16_bits(value)).is_nan());
        } else {
            assert_eq!(f32_to_f16_bits(value), bits, "pattern {bits:#06x}");
        }
    }
}

#[test]
fn values_beyond_the_range_saturate_to_signed_infinity() {
    assert_eq!(f32_to_f16_bits(1.0e9), 0x7c00);
    assert_eq!(f32_to_f16_bits(-1.0e9), 0xfc00);
}

#[test]
fn halfway_cases_round_to_even() {
    // 2049 lies exactly between the binary16 neighbours 2048 and 2050.
    assert_eq!(
        f16_bits_to_f32(f32_to_f16_bits(2049.0)).to_bits(),
        2048.0_f32.to_bits()
    );
    assert_eq!(
        f16_bits_to_f32(f32_to_f16_bits(2051.0)).to_bits(),
        2052.0_f32.to_bits()
    );
}

#[test]
fn a_decimal_with_a_fractional_part_is_not_exact_in_binary16() {
    assert_ne!(
        f16_bits_to_f32(f32_to_f16_bits(31.75 + 0.01)).to_bits(),
        31.76_f32.to_bits()
    );
    assert_eq!(
        f16_bits_to_f32(f32_to_f16_bits(31.75)).to_bits(),
        31.75_f32.to_bits()
    );
}
