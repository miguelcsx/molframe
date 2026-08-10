//! Checked scalar conversions shared by trajectory formats and analyses.

use num_traits::ToPrimitive;

pub(crate) fn f32_from_f64(value: f64) -> Option<f32> {
    value
        .is_finite()
        .then(|| value.to_f32())
        .flatten()
        .filter(|converted| converted.is_finite())
}

pub(crate) fn f32_triplet(values: [f64; 3]) -> Option<[f32; 3]> {
    Some([
        f32_from_f64(values[0])?,
        f32_from_f64(values[1])?,
        f32_from_f64(values[2])?,
    ])
}

pub(crate) fn f64_from_usize(value: usize) -> Option<f64> {
    let converted = value.to_f64()?;
    (converted.to_usize() == Some(value)).then_some(converted)
}

pub(crate) fn f64_from_u64(value: u64) -> Option<f64> {
    let converted = value.to_f64()?;
    (converted.to_u64() == Some(value)).then_some(converted)
}

pub(crate) fn f64_from_i64(value: i64) -> Option<f64> {
    let converted = value.to_f64()?;
    (converted.to_i64() == Some(value)).then_some(converted)
}
