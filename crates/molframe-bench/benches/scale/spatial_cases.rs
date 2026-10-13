//! Spatial-index workloads: neighbour counting, pair enumeration, k-nearest.

use crate::harness::{Axis, Case, Outcome, Run};
use molframe_bench::{Checksum, Tier, uniform_cloud};
use molframe_core::AtomSelection;
use molframe_spatial::{
    KnnOptions, PairQuery, SpatialBackend, SpatialSearchOptions, count_pairs_within, knn_within,
    pairs_within, pairs_within_unsorted,
};

const CUTOFF: f32 = 5.0;
const ALL_TIERS: &[Tier] = &[
    Tier::Tiny,
    Tier::Small,
    Tier::Medium,
    Tier::Large,
    Tier::Huge,
];
const UP_TO_LARGE: &[Tier] = &[Tier::Tiny, Tier::Small, Tier::Medium, Tier::Large];

pub fn cases() -> Vec<Case> {
    vec![
        Case {
            id: "spatial.count_pairs.cell",
            axis: Axis::Atoms,
            tiers: ALL_TIERS,
            run: count_pairs,
        },
        Case {
            id: "spatial.pairs.cell",
            axis: Axis::Atoms,
            tiers: UP_TO_LARGE,
            run: pairs,
        },
        Case {
            id: "spatial.pairs.cell.unsorted",
            axis: Axis::Atoms,
            tiers: UP_TO_LARGE,
            run: pairs_unsorted,
        },
        Case {
            id: "spatial.pairs.kd",
            axis: Axis::Atoms,
            tiers: &[Tier::Tiny, Tier::Small, Tier::Medium],
            run: pairs_kd,
        },
        Case {
            id: "spatial.knn8.kd",
            axis: Axis::Atoms,
            tiers: ALL_TIERS,
            run: knn,
        },
    ]
}

fn selection(count: usize) -> Result<AtomSelection, String> {
    u32::try_from(count)
        .map(AtomSelection::All)
        .map_err(|e| e.to_string())
}

fn count_pairs(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let cloud = uniform_cloud(atoms, 1);
    let all = selection(atoms)?;
    let query = PairQuery {
        positions: &cloud.positions,
        left: &all,
        right: &all,
        cutoff: CUTOFF,
        options: SpatialSearchOptions::with_backend(SpatialBackend::CellList),
        periodic: None,
        context: &run.context,
    };
    let sample = run.measure(
        || count_pairs_within(&query).map_err(|e| e.to_string()),
        |count| *count,
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn enumerate(run: &Run<'_>, atoms: usize, backend: SpatialBackend) -> Result<Outcome, String> {
    enumerate_with(run, atoms, backend, pairs_within)
}

type Enumerator = fn(
    &[[f32; 3]],
    &AtomSelection,
    &AtomSelection,
    f32,
    SpatialBackend,
    Option<&molframe_spatial::PeriodicBox>,
    &molframe_core::ExecutionContext,
)
    -> Result<Vec<molframe_spatial::NeighborPair>, molframe_spatial::SpatialError>;

fn enumerate_with(
    run: &Run<'_>,
    atoms: usize,
    backend: SpatialBackend,
    search: Enumerator,
) -> Result<Outcome, String> {
    let cloud = uniform_cloud(atoms, 1);
    let all = selection(atoms)?;
    let sample = run.measure(
        || {
            search(
                &cloud.positions,
                &all,
                &all,
                CUTOFF,
                backend,
                None,
                &run.context,
            )
            .map_err(|e| e.to_string())
        },
        |pairs| {
            let mut sum = Checksum::default();
            sum.word(pairs.len() as u64);
            for pair in pairs {
                sum.word(u64::from(pair.first) << 32 | u64::from(pair.second));
                sum.word(u64::from(pair.distance_squared.to_bits()));
            }
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}

fn pairs(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    enumerate(run, atoms, SpatialBackend::CellList)
}

fn pairs_unsorted(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    enumerate_with(run, atoms, SpatialBackend::CellList, pairs_within_unsorted)
}

fn pairs_kd(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    enumerate(run, atoms, SpatialBackend::KdTree)
}

fn knn(run: &Run<'_>, atoms: usize) -> Result<Outcome, String> {
    let cloud = uniform_cloud(atoms, 1);
    let all = selection(atoms)?;
    let options = KnnOptions {
        max_distance: None,
        exclude_self: true,
    };
    let sample = run.measure(
        || {
            knn_within(
                &cloud.positions,
                &all,
                &all,
                8,
                options,
                SpatialBackend::KdTree,
                None,
                &run.context,
            )
            .map_err(|e| e.to_string())
        },
        |result| {
            let mut sum = Checksum::default();
            sum.word(result.indices.len() as u64);
            for &index in &result.indices {
                sum.word(u64::from(index));
            }
            sum.value()
        },
    )?;
    Ok(Outcome {
        sample,
        items: atoms as u64,
    })
}
