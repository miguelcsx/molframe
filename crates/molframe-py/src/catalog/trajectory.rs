//! Catalogued operations of the `trajectory` namespace.

use super::{Capability, Domain};
use molframe::Cost;

const D: Domain = Domain("trajectory");

pub(super) const ROWS: &[Capability] = &[
    D.row("read", "path", "Trajectory", "cpu,memory"),
    D.row("rmsd", "positions", "Analysis", "cpu").policy(),
    D.row(
        "pairwise_rmsd",
        "positions",
        "Analysis<array>",
        "cpu,memory",
    )
    .policy(),
    D.row("pca", "positions", "Analysis<Pca>", "cpu,memory")
        .policy(),
    D.row("mean_structure", "positions", "Analysis<array>", "cpu")
        .policy(),
    D.row(
        "diffusion_map",
        "distances",
        "Analysis<DiffusionMap>",
        "cpu,memory",
    )
    .policy(),
    D.row("kmeans", "observations", "Analysis<KMeans>", "cpu")
        .policy(),
    D.row(
        "agglomerative_clustering",
        "distances",
        "Analysis<Clustering>",
        "cpu",
    )
    .policy(),
    D.row("dbscan", "distances", "Analysis<Clustering>", "cpu")
        .policy(),
    D.row("medoid", "distances,members", "Analysis<int>", "cpu")
        .policy(),
    D.row(
        "msd",
        "positions",
        "Analysis<MeanSquaredDisplacement>",
        "cpu",
    )
    .policy(),
    D.row(
        "group_variance",
        "positions,groups",
        "Analysis<list>",
        "cpu",
    )
    .policy(),
    D.row(
        "block_convergence",
        "values",
        "Analysis<Convergence>",
        "cpu",
    )
    .policy(),
    D.row(
        "harmonic_similarity",
        "first,second",
        "Analysis<HarmonicSimilarity>",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "population_similarity",
        "first,second",
        "Analysis<float>",
        "cpu",
    )
    .policy(),
    D.row(
        "path_similarity",
        "first,second",
        "PathSimilarity",
        "cpu,memory",
    ),
    D.row("read_chunks", "path", "TrajectoryChunks", "cpu,memory")
        .cost(Cost::Decode),
    D.row("rmsd_stream", "path,reference", "FrameSeries", "cpu,memory")
        .cost(Cost::Decode),
    D.row("rmsf_stream", "path", "array", "cpu,memory")
        .cost(Cost::Decode),
    D.row(
        "sasa_stream",
        "path,topology,probe,n_points,hetero,metals",
        "Analysis<FrameSeries>",
        "cpu,memory",
    )
    .cost(Cost::Decode)
    .policy(),
    D.row(
        "contact_counts_stream",
        "path,cutoff",
        "FrameSeries",
        "cpu,memory",
    )
    .cost(Cost::Decode),
    D.row("run_analysis_stream", "path,callback", "int", "cpu,memory")
        .cost(Cost::Decode),
];
