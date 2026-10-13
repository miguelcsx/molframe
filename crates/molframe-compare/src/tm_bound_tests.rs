use super::TmBound;
use crate::tm_score;

/// A deterministic wandering chain.
pub(crate) fn chain(seed: u32, atoms: usize) -> Vec<[f32; 3]> {
    let mut state = seed.wrapping_mul(2_654_435_761).wrapping_add(12_345);
    let mut next = move || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        f32::from(u16::try_from(state >> 16).unwrap_or(0)) / 65_536.0 - 0.5
    };
    let mut point = [0.0_f32; 3];
    (0..atoms)
        .map(|_| {
            point = [
                point[0] + 3.0 * next() + 1.5,
                point[1] + 3.0 * next(),
                point[2] + 3.0 * next(),
            ];
            point
        })
        .collect()
}

/// `chain` with every point nudged by up to `jitter`.
pub(crate) fn nudged(base: &[[f32; 3]], seed: u32, jitter: f32) -> Vec<[f32; 3]> {
    let mut state = seed.wrapping_mul(747_796_405).wrapping_add(2_891_336_453);
    let mut next = move || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        (f32::from(u16::try_from(state >> 16).unwrap_or(0)) / 65_536.0 - 0.5) * 2.0 * jitter
    };
    base.iter()
        .map(|p| [p[0] + next(), p[1] + next(), p[2] + next()])
        .collect()
}

#[test]
fn bound_never_falls_below_the_score() {
    let mut bound = TmBound::default();
    for (length, seed) in [(5, 1), (12, 2), (40, 3), (120, 4), (250, 5), (600, 6)] {
        let base = chain(seed, length);
        let others = [
            nudged(&base, seed, 0.3),
            nudged(&base, seed + 9, 2.0),
            nudged(&base, seed + 17, 6.0),
            chain(seed + 100, length),
        ];
        for other in &others {
            let score = tm_score(other, &base).expect("scores");
            let upper = bound.bound(other, &base);
            assert!(
                upper >= score,
                "bound {upper} below score {score} at length {length}"
            );
            assert!(upper <= 1.0);
        }
    }
}

#[test]
fn identical_structures_are_not_bounded_below_one() {
    let base = chain(7, 80);
    assert!(TmBound::default().bound(&base, &base) >= 1.0 - 1e-9);
}

#[test]
fn unrelated_chains_are_bounded_well_under_one() {
    let a = chain(11, 200);
    let b = chain(12, 200);
    let upper = TmBound::default().bound(&a, &b);
    assert!(
        upper < 0.95,
        "bound {upper} should separate unrelated chains"
    );
}

#[test]
fn bound_is_stable_across_scratch_reuse() {
    let a = chain(21, 90);
    let b = nudged(&a, 3, 1.0);
    let mut reused = TmBound::default();
    let _ = reused.bound(&chain(1, 300), &chain(2, 300));
    assert_eq!(
        reused.bound(&a, &b).to_bits(),
        TmBound::default().bound(&a, &b).to_bits()
    );
}
