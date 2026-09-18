/// Density-based clustering with explicit neighbourhood and core-size policy.
///
/// # Errors
///
/// Returns [`EnsembleGeometryError`] for an invalid matrix, epsilon, or core size.
pub fn dbscan_clustering(
    distances: &EnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
) -> Result<Clustering, EnsembleGeometryError> {
    validate_matrix(distances)?;

    if !epsilon.is_finite() || epsilon < 0.0 || minimum_points == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }

    let size = distances.size;
    let mut labels = vec![None; size];
    let mut visited = vec![false; size];
    let mut queued_epoch = vec![usize::MAX; size];

    let mut neighbours = reusable_index_buffer(size)?;
    let mut queue = reusable_index_buffer(size)?;
    let mut scratch = reusable_index_buffer(size)?;

    let mut cluster_count = 0usize;

    for point in 0..size {
        let visited_slot = visited
            .get_mut(point)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        if *visited_slot {
            continue;
        }

        *visited_slot = true;
        collect_neighbourhood(distances, point, epsilon, &mut neighbours)?;

        if neighbours.len() < minimum_points {
            continue;
        }

        let mut expansion = ClusterExpansion {
            distances,
            epsilon,
            minimum_points,
            label: cluster_count,
            visited: &mut visited,
            labels: &mut labels,
            queued_epoch: &mut queued_epoch,
        };

        expansion.run(&neighbours, &mut queue, &mut scratch)?;

        cluster_count = cluster_count
            .checked_add(1)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
    }

    let mut members = vec![Vec::new(); cluster_count];

    for (point, label) in labels.iter().copied().enumerate() {
        let Some(label) = label else {
            continue;
        };

        let cluster = members
            .get_mut(label)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        cluster.push(point);
    }

    let medoids = cluster_medoids(distances, &members)?;

    Ok(Clustering {
        labels,
        members,
        medoids,
    })
}
