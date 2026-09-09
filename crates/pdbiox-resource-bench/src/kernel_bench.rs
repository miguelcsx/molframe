//! Deterministic in-process kernel probes: trajectory, pore, and alignment.

use super::{ResourceRecord, measure_case};
use std::hint::black_box;

pub(super) fn run_trajectory_contacts() -> Result<ResourceRecord, String> {
    let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Tiny);
    let trajectory = pdbiox::traj::Trajectory::from_frames(
        (0_u16..64)
            .map(|frame| pdbiox::traj::Frame {
                positions: structure
                    .positions()
                    .iter()
                    .map(|&[x, y, z]| [x + f32::from(frame) * 1.0e-4, y, z])
                    .collect(),
            })
            .collect(),
    )
    .map_err(|error| format!("trajectory construction failed: {error}"))?;
    let policy = pdbiox::AnalysisPolicy::default();
    let kernel = pdbiox::analysis::contacts_kernel(3.0, pdbiox::SpatialBackend::Auto);
    measure_case("trajectory_contacts", || {
        let analysis = pdbiox::analysis::analyse_trajectory(
            &structure,
            &trajectory,
            &policy,
            &kernel,
            &pdbiox::ExecutionContext::builder()
                .worker_budget(4)
                .build()
                .map_err(|error| format!("execution context failed: {error}"))?,
        )
        .map_err(|error| format!("trajectory analysis failed: {error}"))?;
        black_box(analysis.value);
        Ok(64)
    })
}

pub(super) fn run_pore_profile_100000() -> Result<ResourceRecord, String> {
    let mut positions = Vec::with_capacity(100_000);
    for ring in 0_u16..1_000 {
        let z = -50.0 + 100.0 * f32::from(ring) / 999.0;
        for atom in 0_u16..100 {
            let angle = std::f32::consts::TAU * f32::from(atom) / 100.0;
            let radius = 9.0 + 0.5 * (angle * 5.0 + z * 0.1).sin();
            positions.push([radius * angle.cos(), radius * angle.sin(), z]);
        }
    }
    let radii = vec![1.7; positions.len()];
    let options = pdbiox::analysis::PoreProfileOptions::new(
        ([0.0; 3], [0.0, 0.0, 1.0]),
        -10.0,
        10.0,
        16,
        5.0,
        1.0,
        1.4,
    );
    measure_case("pore_profile_100000", || {
        let profile = pdbiox::analysis::pore_profile(&positions, &radii, options)
            .map_err(|error| format!("pore profile failed: {error}"))?;
        Ok(profile.iter().fold(0_u64, |digest, sample| {
            digest.wrapping_add(u64::from(sample.radius.to_bits()))
        }))
    })
}

pub(super) fn run_msa_eight_512() -> Result<ResourceRecord, String> {
    let sequences: [Vec<u8>; 8] = std::array::from_fn(|sequence_index| {
        (0..512 + sequence_index * 3)
            .map(|position| {
                const ALPHABET: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
                let base = ALPHABET[(position * 17 + position / 7) % ALPHABET.len()];
                if position >= sequence_index + 11
                    && (position - sequence_index - 11).is_multiple_of(47 + sequence_index)
                {
                    b'Y'
                } else {
                    base
                }
            })
            .collect()
    });
    let views = sequences.each_ref().map(Vec::as_slice);
    measure_case("msa_eight_512", || {
        let alignment = pdbiox::seq::progressive_msa(
            &views,
            pdbiox::seq::MsaOptions::progressive(pdbiox::seq::Scoring::simple()),
        )
        .map_err(|error| format!("MSA failed: {error}"))?;
        alignment_digest(&alignment)
    })
}

pub(super) fn run_msa_sixty_four_128() -> Result<ResourceRecord, String> {
    let sequences = (0..64)
        .map(|sequence_index: usize| {
            (0..128)
                .map(|position: usize| {
                    const ALPHABET: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
                    if position >= sequence_index % 23
                        && (position - sequence_index % 23).is_multiple_of(29 + sequence_index % 7)
                    {
                        b'W'
                    } else {
                        ALPHABET[(position * 17 + position / 7) % ALPHABET.len()]
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let views = sequences.iter().map(Vec::as_slice).collect::<Vec<_>>();
    measure_case("msa_sixty_four_128", || {
        let alignment = pdbiox::seq::progressive_msa(
            &views,
            pdbiox::seq::MsaOptions::progressive(pdbiox::seq::Scoring::simple()),
        )
        .map_err(|error| format!("MSA failed: {error}"))?;
        alignment_digest(&alignment)
    })
}

fn alignment_digest(alignment: &[Vec<u8>]) -> Result<u64, String> {
    alignment.iter().try_fold(0_u64, |total, row| {
        u64::try_from(row.len())
            .ok()
            .and_then(|length| total.checked_add(length))
            .ok_or_else(|| "MSA digest exceeds u64".to_owned())
    })
}
