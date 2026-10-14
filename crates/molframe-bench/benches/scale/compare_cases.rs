//! Structure-comparison workloads: fitted scores, lDDT, all-vs-all.

use crate::harness::{Axis, Case, Outcome, Run};
use molframe_bench::{Checksum, Tier, chain, perturbed, rotated};
use molframe_compare::{
    TmPrefilter, distance_difference, gdt_ts, tm_score, tm_score_all_pairs, tm_screen_all_pairs,
};

const ALL_TIERS: &[Tier] = &[
    Tier::Tiny,
    Tier::Small,
    Tier::Medium,
    Tier::Large,
    Tier::Huge,
];
const RESIDUES: usize = 250;

pub fn cases() -> Vec<Case> {
    vec![
        Case {
            id: "compare.tm_score",
            axis: Axis::Atoms,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium, Tier::Large],
            run: tm_case,
        },
        Case {
            id: "compare.gdt_ts",
            axis: Axis::Atoms,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium, Tier::Large],
            run: gdt_case,
        },
        Case {
            id: "compare.distance_difference",
            axis: Axis::Atoms,
            tiers: ALL_TIERS,
            run: distance_difference_case,
        },
        Case {
            id: "compare.all_vs_all_tm.250_res",
            axis: Axis::Structures,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium],
            run: all_vs_all,
        },
        Case {
            id: "compare.all_pairs_tm.250_res",
            axis: Axis::Structures,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium],
            run: all_pairs_batch,
        },
        Case {
            id: "compare.mixed_all_pairs_tm.250_res",
            axis: Axis::Structures,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium],
            run: mixed_unfiltered,
        },
        Case {
            id: "compare.mixed_screen_exact_tm.250_res",
            axis: Axis::Structures,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium],
            run: mixed_exact,
        },
        Case {
            id: "compare.mixed_screen_coarse_tm.250_res",
            axis: Axis::Structures,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium],
            run: mixed_coarse,
        },
    ]
}

fn pair(atoms: usize) -> (Vec<[f32; 3]>, Vec<[f32; 3]>) {
    let reference = chain(atoms, 3);
    let model = perturbed(&rotated(&reference, 4), 0.2);
    (reference, model)
}

fn tm_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let (reference, model) = pair(atoms);
    let sample = run.measure(
        || tm_score(&model, &reference).map_err(|e| e.to_string()),
        |value| crate::harness::real_fingerprint(*value),
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn gdt_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let (reference, model) = pair(atoms);
    let sample = run.measure(
        || gdt_ts(&model, &reference).map_err(|e| e.to_string()),
        |value| crate::harness::real_fingerprint(*value),
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn distance_difference_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let (reference, model) = pair(atoms);
    let sample = run.measure(
        || distance_difference(&model, &reference, 15.0, &run.context).map_err(|e| e.to_string()),
        |value| crate::harness::real_fingerprint(*value),
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

/// Every unordered pair of `count` 250-residue models through [`tm_score`],
/// one pair per serial iteration: the cost a caller pays without a batch API.
fn all_vs_all(run: &Run<'_>, count: usize) -> Result<Outcome, String> {
    let models = models(count);
    let pairs = (count * (count - 1) / 2) as u64;
    let sample = run.measure(
        || {
            let mut values = Vec::with_capacity(pairs as usize);
            for i in 0..count {
                for j in (i + 1)..count {
                    values.push(tm_score(&models[i], &models[j]).map_err(|e| e.to_string())?);
                }
            }
            Ok(values)
        },
        |values: &Vec<f64>| {
            let mut sum = Checksum::default();
            sum.reals(values);
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: pairs,
    })
}

fn models(count: usize) -> Vec<Vec<[f32; 3]>> {
    let base = chain(RESIDUES, 3);
    (0..count)
        .map(|index| perturbed(&rotated(&base, 100 + index as u64), 0.5))
        .collect()
}

/// The same pairs through the batch API, parallel across pairs.
fn all_pairs_batch(run: &Run<'_>, count: usize) -> Result<Outcome, String> {
    let owned = models(count);
    let views: Vec<&[[f32; 3]]> = owned.iter().map(Vec::as_slice).collect();
    let pairs = (count * (count - 1) / 2) as u64;
    let sample = run.measure(
        || tm_score_all_pairs(&views, &run.context).map_err(|e| e.to_string()),
        |values: &Vec<f64>| {
            let mut sum = Checksum::default();
            sum.reals(values);
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: pairs,
    })
}

/// Near-copies in families of four, so most pairs are unrelated folds.
fn mixed_models(count: usize) -> Vec<Vec<[f32; 3]>> {
    (0..count)
        .map(|index| {
            let family = (index / 4) as u64;
            let base = chain(RESIDUES, 1000 + family);
            perturbed(&rotated(&base, 100 + index as u64), 0.5)
        })
        .collect()
}

/// The exhaustive batch over the mixed library: the reference to screen.
fn mixed_unfiltered(run: &Run<'_>, count: usize) -> Result<Outcome, String> {
    let owned = mixed_models(count);
    let views: Vec<&[[f32; 3]]> = owned.iter().map(Vec::as_slice).collect();
    let pairs = (count * (count - 1) / 2) as u64;
    let sample = run.measure(
        || tm_score_all_pairs(&views, &run.context).map_err(|e| e.to_string()),
        |values: &Vec<f64>| {
            let mut sum = Checksum::default();
            sum.reals(values);
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: pairs,
    })
}

fn mixed_exact(run: &Run<'_>, count: usize) -> Result<Outcome, String> {
    let filter = TmPrefilter::new(0.5).map_err(|e| e.to_string())?;
    mixed_screen(run, count, &filter)
}

fn mixed_coarse(run: &Run<'_>, count: usize) -> Result<Outcome, String> {
    let filter = TmPrefilter::new(0.5)
        .and_then(|f| f.with_coarse_ratio(0.5))
        .map_err(|e| e.to_string())?;
    mixed_screen(run, count, &filter)
}

/// The mixed library through the screen; items are the pairs considered.
fn mixed_screen(run: &Run<'_>, count: usize, filter: &TmPrefilter) -> Result<Outcome, String> {
    let owned = mixed_models(count);
    let views: Vec<&[[f32; 3]]> = owned.iter().map(Vec::as_slice).collect();
    let pairs = (count * (count - 1) / 2) as u64;
    if let Ok(screen) = tm_screen_all_pairs(&views, filter, &run.context) {
        eprintln!(
            "scale: screen n={count} evaluated={} bound={} coarse={}",
            screen.evaluated, screen.pruned_bound, screen.pruned_coarse
        );
    }
    let sample = run.measure(
        || tm_screen_all_pairs(&views, filter, &run.context).map_err(|e| e.to_string()),
        |screen| {
            let mut sum = Checksum::default();
            sum.reals(&screen.scores);
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: pairs,
    })
}
