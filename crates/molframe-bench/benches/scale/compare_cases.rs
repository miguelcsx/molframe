//! Structure-comparison workloads: fitted scores, lDDT, all-vs-all.

use crate::harness::{Axis, Case, Outcome, Run};
use molframe_bench::{Checksum, Tier, chain, perturbed, rotated};
use molframe_compare::{gdt_ts, lddt, tm_score};

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
            id: "compare.lddt",
            axis: Axis::Atoms,
            tiers: ALL_TIERS,
            run: lddt_case,
        },
        Case {
            id: "compare.all_vs_all_tm.250_res",
            axis: Axis::Structures,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium],
            run: all_vs_all,
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
        crate::harness::real_fingerprint,
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
        crate::harness::real_fingerprint,
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn lddt_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let (reference, model) = pair(atoms);
    let sample = run.measure(
        || lddt(&model, &reference, 15.0, &run.context).map_err(|e| e.to_string()),
        crate::harness::real_fingerprint,
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

/// Every unordered pair of `count` 250-residue models through [`tm_score`],
/// one pair per serial iteration: the cost a caller pays without a batch API.
fn all_vs_all(run: &Run<'_>, count: usize) -> Result<Outcome, String> {
    let base = chain(RESIDUES, 3);
    let models: Vec<Vec<[f32; 3]>> = (0..count)
        .map(|index| perturbed(&rotated(&base, 100 + index as u64), 0.5))
        .collect();
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
