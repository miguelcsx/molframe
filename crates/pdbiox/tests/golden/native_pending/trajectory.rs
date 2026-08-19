use pdbiox::{AnalysisPolicy, AtomSelection, SpatialBackend};

#[test]
fn gw_024_reuses_a_neighbour_list_across_trajectory_frames() {
    let positions = [
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.1, 0.0, 0.0], [1.1, 0.0, 0.0]],
    ];
    let selection = AtomSelection::All(2);
    let options = pdbiox::analysis::RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 2.0,
        bins: 4,
        volume: 100.0,
        backend: SpatialBackend::BruteForce,
    };
    let mut neighbours = pdbiox::traj::FrameNeighborList::new([0_u32, 1], [0_u32, 1], 2.0, 0.5);
    let mut radial_counts = Vec::new();
    for frame in positions {
        let _ = neighbours
            .pairs_positions(&frame, None)
            .unwrap_or_else(|error| panic!("neighbour-list query failed: {error}"));
        radial_counts.push(
            pdbiox::analysis::radial_distribution(&frame, &selection, &selection, options, None)
                .unwrap_or_else(|error| panic!("RDF failed: {error}"))
                .iter()
                .map(|bin| bin.count)
                .sum::<u64>(),
        );
    }
    assert_eq!(radial_counts, [1, 1]);
    assert_eq!(neighbours.statistics().frames, 2);
    assert_eq!(neighbours.statistics().rebuilds, 1);
}

#[test]
fn gw_025_extracts_selected_frames_and_writes_them_losslessly() {
    let source = r"2
frame-0
C 0 0 0
N 1 0 0
2
frame-1
C 1 0 0
N 2 0 0
2
frame-2
C 2 0 0
N 3 0 0
";
    let frames = pdbiox::traj::parse_xyz(source).unwrap_or_else(|| panic!("XYZ parse failed"));
    let selected: Vec<_> = [0_usize, 2]
        .iter()
        .map(|index| frames[*index].clone())
        .collect();
    let rendered = pdbiox::traj::write_xyz(&selected);
    let round_trip =
        pdbiox::traj::parse_xyz(&rendered).unwrap_or_else(|| panic!("selected XYZ parse failed"));
    assert_eq!(round_trip, selected);
    assert_eq!(round_trip.len(), 2);
    assert_eq!(round_trip[1].comment, "frame-2");
}

#[test]
fn gw_026_parallel_trajectory_analysis_is_byte_identical() {
    let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Tiny);
    let trajectory = pdbiox::traj::Trajectory::from_frames(
        (0_u16..1_025)
            .map(|frame| pdbiox::traj::Frame {
                positions: structure
                    .positions()
                    .iter()
                    .map(|&[x, y, z]| [x + f32::from(frame) * 1.0e-4, y, z])
                    .collect(),
            })
            .collect(),
    );
    let policy = AnalysisPolicy::default();
    let kernel = pdbiox::analysis::contacts_kernel(3.0, SpatialBackend::Auto);
    let results: Vec<_> = [1_usize, 2, 4, 16]
        .into_iter()
        .map(|workers| {
            pdbiox::analysis::analyse_trajectory(&structure, &trajectory, &policy, &kernel, workers)
                .unwrap_or_else(|error| panic!("parallel analysis failed at {workers}: {error}"))
        })
        .collect();
    assert!(results.windows(2).all(|pair| {
        pair[0].value == pair[1].value
            && pair[0].coverage == pair[1].coverage
            && pair[0].provenance.fingerprint() == pair[1].provenance.fingerprint()
    }));
}
