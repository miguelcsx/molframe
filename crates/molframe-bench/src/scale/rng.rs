//! A tiny seeded generator, so inputs are bit-identical on every machine.

/// `SplitMix64`: one `u64` of state, full period, no dependencies.
#[derive(Clone, Debug)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// A generator whose whole sequence is a function of `seed`.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform value in `[0, 1)` with 24 bits of resolution.
    pub fn unit(&mut self) -> f32 {
        let bits = molframe_core::numeric::u64_to_u32(self.next_u64() >> 40);
        let numerator = f32::from(molframe_core::numeric::low_u16(bits >> 8)) * 256.0
            + f32::from(molframe_core::numeric::low_u16(bits & 0xFF));
        numerator / 16_777_216.0
    }

    /// An approximately normal value (Irwin-Hall, twelve uniforms), mean 0, sd 1.
    pub fn normal(&mut self) -> f32 {
        let mut sum = 0.0;
        for _ in 0..12 {
            sum += self.unit();
        }
        sum - 6.0
    }

    /// A uniformly distributed unit vector.
    pub fn direction(&mut self) -> [f32; 3] {
        loop {
            let v = [
                2.0 * self.unit() - 1.0,
                2.0 * self.unit() - 1.0,
                2.0 * self.unit() - 1.0,
            ];
            let n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            if n2 > 1.0e-4 && n2 <= 1.0 {
                let n = n2.sqrt();
                return [v[0] / n, v[1] / n, v[2] / n];
            }
        }
    }
}
