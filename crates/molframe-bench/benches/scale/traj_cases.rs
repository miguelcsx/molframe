//! Trajectory workloads: streaming metrics, PCA, per-frame contacts.

use crate::harness::{Axis, Case, Outcome, Run};
use molframe_bench::{Checksum, Tier, ensemble};
use molframe_spatial::{SpatialBackend, SpatialSearchOptions};
use molframe_traj::{
    CartesianFit, Frame, FrameAlignment, MemoryReader, Trajectory, cartesian_pca_view,
    contact_counts_stream, rmsd_stream, rmsf_stream,
};

const ATOMS: usize = 1_000;
const FRAME_BYTES: usize = ATOMS * 12 * 4;
const FRAME_TIERS: &[Tier] = &[
    Tier::Tiny,
    Tier::Small,
    Tier::Medium,
    Tier::Large,
    Tier::Huge,
];

pub fn cases() -> Vec<Case> {
    vec![
        Case {
            id: "traj.rmsd_stream.rigid.1k_atoms",
            axis: Axis::Frames,
            tiers: FRAME_TIERS,
            run: rmsd_case,
        },
        Case {
            id: "traj.rmsf_stream.1k_atoms",
            axis: Axis::Frames,
            tiers: FRAME_TIERS,
            run: rmsf_case,
        },
        Case {
            id: "traj.contacts_stream.1k_atoms",
            axis: Axis::Frames,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium, Tier::Large],
            run: contacts_case,
        },
        Case {
            id: "traj.cartesian_pca.100_atoms",
            axis: Axis::Frames,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium, Tier::Large],
            run: pca_case,
        },
    ]
}

fn trajectory(atoms: usize, frames: usize) -> Result<Trajectory, String> {
    let ens = ensemble(atoms, frames, 0.8, 5);
    let frames = (0..frames)
        .map(|f| Frame {
            positions: ens.frame(f).to_vec(),
        })
        .collect();
    Trajectory::from_frames(frames).map_err(|e| e.to_string())
}

fn rmsd_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let trajectory = trajectory(ATOMS, frames)?;
    let reference = trajectory.positions()[..ATOMS].to_vec();
    let sample = run.measure(
        || {
            let mut reader = MemoryReader::new(&trajectory);
            let mut values = Vec::with_capacity(frames);
            rmsd_stream(
                &mut reader,
                &reference,
                FrameAlignment::Rigid,
                &run.context,
                FRAME_BYTES,
                |_, _, value| {
                    values.push(value);
                    Ok(())
                },
            )
            .map_err(|e| e.to_string())?;
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
        items: frames as u64,
    })
}

fn rmsf_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let trajectory = trajectory(ATOMS, frames)?;
    let sample = run.measure(
        || {
            let mut reader = MemoryReader::new(&trajectory);
            rmsf_stream(&mut reader, &run.context, FRAME_BYTES).map_err(|e| e.to_string())
        },
        |values| {
            let mut sum = Checksum::default();
            sum.reals(values);
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: (frames * ATOMS) as u64,
    })
}

fn contacts_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let trajectory = trajectory(ATOMS, frames)?;
    let options = SpatialSearchOptions::with_backend(SpatialBackend::CellList);
    let sample = run.measure(
        || {
            let mut reader = MemoryReader::new(&trajectory);
            let mut counts = Vec::with_capacity(frames);
            contact_counts_stream(
                &mut reader,
                6.0,
                options,
                &run.context,
                FRAME_BYTES,
                |_, _, count| {
                    counts.push(count);
                    Ok(())
                },
            )
            .map_err(|e| e.to_string())?;
            Ok(counts)
        },
        |counts: &Vec<u64>| {
            let mut sum = Checksum::default();
            for &count in counts {
                sum.word(count);
            }
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: frames as u64,
    })
}

fn pca_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let atoms = 100;
    let ens = ensemble(atoms, frames, 0.8, 5);
    let view =
        molframe_traj::FrameView::new(&ens.positions, frames, atoms).map_err(|e| e.to_string())?;
    let sample = run.measure(
        || {
            cartesian_pca_view(view, CartesianFit::None, 3, usize::MAX / 4)
                .map_err(|e| format!("{e:?}"))
        },
        |pca| {
            let mut sum = Checksum::default();
            sum.reals(&pca.eigenvalues);
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: frames as u64,
    })
}
