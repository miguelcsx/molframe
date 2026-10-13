//! Coordinate-kernel workloads: fitting, fluctuation, moments, distances.

use crate::harness::{Axis, Case, Outcome, Run};
use molframe_bench::{Checksum, Tier, chain, ensemble, perturbed, rotated};
use molframe_geom::{
    centroid, condensed_distances_f32_with_context, distance_matrix_with_context, inertia_tensor,
    rmsd, rmsd_after_fit, rmsf, superpose,
};

const ALL_TIERS: &[Tier] = &[
    Tier::Tiny,
    Tier::Small,
    Tier::Medium,
    Tier::Large,
    Tier::Huge,
];

pub fn cases() -> Vec<Case> {
    vec![
        Case {
            id: "geom.rmsd",
            axis: Axis::Atoms,
            tiers: ALL_TIERS,
            run: rmsd_case,
        },
        Case {
            id: "geom.superpose",
            axis: Axis::Atoms,
            tiers: ALL_TIERS,
            run: superpose_case,
        },
        Case {
            id: "geom.inertia_tensor",
            axis: Axis::Atoms,
            tiers: ALL_TIERS,
            run: inertia_case,
        },
        Case {
            id: "geom.rmsf.1k_atoms",
            axis: Axis::Frames,
            tiers: ALL_TIERS,
            run: rmsf_case,
        },
        Case {
            id: "geom.distance_matrix",
            axis: Axis::Atoms,
            tiers: &[Tier::Tiny, Tier::Small],
            run: matrix_case,
        },
        Case {
            id: "geom.condensed_f32",
            axis: Axis::Atoms,
            tiers: &[Tier::Tiny, Tier::Small],
            run: condensed_case,
        },
    ]
}

fn pair(atoms: usize) -> (Vec<[f32; 3]>, Vec<[f32; 3]>) {
    let reference = chain(atoms, 3);
    let model = perturbed(&rotated(&reference, 4), 0.2);
    (reference, model)
}

fn rmsd_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let (reference, model) = pair(atoms);
    let sample = run.measure(
        || rmsd(&model, &reference).map_err(|e| format!("{e:?}")),
        crate::harness::real_fingerprint,
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn superpose_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let (reference, model) = pair(atoms);
    let sample = run.measure(
        || rmsd_after_fit(&model, &reference).map_err(|e| format!("{e:?}")),
        crate::harness::real_fingerprint,
    )?;
    let _ = superpose;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn inertia_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let reference = chain(atoms, 3);
    let masses = vec![12.0; atoms];
    let sample = run.measure(
        || {
            let centre = centroid(&reference).ok_or("empty")?;
            let tensor = inertia_tensor(&reference, &masses).ok_or("empty")?;
            Ok((centre, tensor))
        },
        |(centre, tensor): &([f64; 3], [[f64; 3]; 3])| {
            let mut sum = Checksum::default();
            sum.reals(centre);
            for row in tensor {
                sum.reals(row);
            }
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn rmsf_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let atoms = 1_000;
    let ens = ensemble(atoms, frames, 0.8, 5);
    let views: Vec<&[[f32; 3]]> = (0..frames).map(|f| ens.frame(f)).collect();
    let sample = run.measure(
        || rmsf(&views).map_err(|e| format!("{e:?}")),
        |values| {
            let mut sum = Checksum::default();
            sum.reals(values);
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: (frames * atoms) as u64,
    })
}

fn matrix_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let points = chain(atoms, 3);
    let sample = run.measure(
        || distance_matrix_with_context(&points, &run.context).map_err(|e| format!("{e:?}")),
        |matrix| {
            let mut sum = Checksum::default();
            sum.reals(matrix.as_slice());
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: (atoms * atoms) as u64,
    })
}

fn condensed_case(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let points = chain(atoms, 3);
    let sample = run.measure(
        || {
            condensed_distances_f32_with_context(&points, &run.context)
                .map_err(|e| format!("{e:?}"))
        },
        |values| {
            let mut sum = Checksum::default();
            for value in values {
                sum.word(u64::from(value.to_bits()));
            }
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: (atoms * atoms / 2) as u64,
    })
}
