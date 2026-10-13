use super::{PaeKind, PaePolicy, PaeStorage};
use crate::ConfidenceError;
use crate::PrecisionProbe;

fn probe(values: &[f32]) -> PrecisionProbe {
    let mut probe = PrecisionProbe::new();
    for value in values {
        probe.observe(*value);
    }
    probe
}

#[test]
fn whole_number_error_maps_to_one_byte_cells() {
    let layout = PaePolicy::new()
        .layout(&probe(&[0.0, 12.0, 31.0]), 4)
        .expect("layout");
    assert_eq!((layout.kind, layout.denominator), (PaeKind::U8, 1));
}

#[test]
fn quarter_steps_beyond_the_byte_range_fit_in_sixteen_bits() {
    let layout = PaePolicy::new()
        .layout(&probe(&[0.25, 31.75]), 4)
        .expect("layout");
    assert_eq!((layout.kind, layout.denominator), (PaeKind::U16, 100));
    let wide = PaePolicy::new()
        .layout(&probe(&[0.25, 3100.75]), 4)
        .expect("layout");
    assert_eq!(wide.kind, PaeKind::F32);
}

#[test]
fn values_without_short_decimals_fall_back_to_binary32() {
    let layout = PaePolicy::new()
        .layout(&probe(&[0.123_456_7, 1.0e-7]), 4)
        .expect("layout");
    assert_eq!(layout.kind, PaeKind::F32);
}

#[test]
fn forcing_a_narrow_layout_chooses_a_coarser_denominator() {
    let policy = PaePolicy::new().with_storage(PaeStorage::Quantised(PaeKind::U8));
    let layout = policy.layout(&probe(&[0.123_4, 31.75]), 4).expect("layout");
    assert_eq!(layout.kind, PaeKind::U8);
    assert_eq!(layout.denominator, 8);
}

#[test]
fn a_matrix_beyond_the_budget_is_refused_before_allocation() {
    let policy = PaePolicy::new().with_max_bytes(100);
    let found = policy.layout(&probe(&[1.0]), 11);
    assert_eq!(
        found,
        Err(ConfidenceError::Budget {
            required: 121,
            limit: 100
        })
    );
}
