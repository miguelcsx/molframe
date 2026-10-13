//! Scaling-exponent fitting and order-sensitive result fingerprints.

/// A least-squares line through `(ln size, ln time)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fit {
    /// The slope: time grows as `size.powf(exponent)`.
    pub exponent: f64,
    /// The coefficient of determination of the log-log line.
    pub r_squared: f64,
}

/// Fits `time = c * size^k` to positive `(size, time)` samples.
///
/// Returns `None` with fewer than two distinct sizes or any non-positive value.
#[must_use]
pub fn fit_exponent(samples: &[(f64, f64)]) -> Option<Fit> {
    if samples.len() < 2 || samples.iter().any(|&(x, y)| x <= 0.0 || y <= 0.0) {
        return None;
    }
    let points: Vec<(f64, f64)> = samples.iter().map(|&(x, y)| (x.ln(), y.ln())).collect();
    let count = f64::from(u32::try_from(points.len()).ok()?);
    let mean_x = points.iter().map(|p| p.0).sum::<f64>() / count;
    let mean_y = points.iter().map(|p| p.1).sum::<f64>() / count;
    let sxx = points.iter().map(|p| (p.0 - mean_x).powi(2)).sum::<f64>();
    let syy = points.iter().map(|p| (p.1 - mean_y).powi(2)).sum::<f64>();
    let sxy = points
        .iter()
        .map(|p| (p.0 - mean_x) * (p.1 - mean_y))
        .sum::<f64>();
    if sxx <= 0.0 {
        return None;
    }
    let r_squared = if syy <= 0.0 {
        1.0
    } else {
        sxy * sxy / (sxx * syy)
    };
    Some(Fit {
        exponent: sxy / sxx,
        r_squared,
    })
}

/// An order-sensitive 64-bit FNV-1a fingerprint of a result, bit for bit.
///
/// Two runs that disagree in a single bit of any folded value, or in the order
/// of values, disagree here; that is the determinism contract across worker
/// counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Checksum(u64);

impl Default for Checksum {
    fn default() -> Self {
        Self(0xCBF2_9CE4_8422_2325)
    }
}

impl Checksum {
    /// Folds raw bytes.
    pub fn bytes(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 = (self.0 ^ u64::from(byte)).wrapping_mul(0x0100_0000_01B3);
        }
    }

    /// Folds an integer.
    pub fn word(&mut self, value: u64) {
        self.bytes(&value.to_le_bytes());
    }

    /// Folds the exact bits of a double.
    pub fn real(&mut self, value: f64) {
        self.word(value.to_bits());
    }

    /// Folds the exact bits of every double in a slice.
    pub fn reals(&mut self, values: &[f64]) {
        for &value in values {
            self.real(value);
        }
    }

    /// The fingerprint.
    #[must_use]
    pub const fn value(self) -> u64 {
        self.0
    }
}
