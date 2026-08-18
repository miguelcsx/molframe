/// Agglomerative clustering to an explicit number of clusters.
///
/// Ties are resolved by the lexicographically smallest pair of current member
/// lists, making results independent of hash iteration or thread scheduling.
///
/// # Errors
///
/// Returns [`EnsembleGeometryError`] for an invalid matrix or cluster count.
pub fn agglomerative_clustering(
    distances: &EnsembleDistanceMatrix,
    cluster_count: usize,
    linkage: Linkage,
) -> Result<Clustering, EnsembleGeometryError> {
    validate_matrix(distances)?;

    if cluster_count == 0 || cluster_count > distances.size {
        return Err(EnsembleGeometryError::InvalidParameter);
    }

    let size = distances.size;
    let mut members: Vec<Vec<usize>> = (0..size).map(|index| vec![index]).collect();

    if cluster_count == size {
        return make_clustering(distances, members, &[]);
    }

    let mut pair_distances = PairDistances::from_matrix(distances)?;
    let mut heap = initial_merge_heap(&pair_distances)?;

    let mut active = vec![true; size];
    let mut versions = vec![0usize; size];
    let mut sizes = vec![1usize; size];
    let mut remaining = size;

    while remaining > cluster_count {
        let candidate = pop_current_candidate(&mut heap, &active, &versions)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        let left = candidate.left;
        let right = candidate.right;

        let left_size = sizes
            .get(left)
            .copied()
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        let right_size = sizes
            .get(right)
            .copied()
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        let merged_size = left_size
            .checked_add(right_size)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        let next_version = versions
            .get(left)
            .copied()
            .and_then(|version| version.checked_add(1))
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        merge_member_slots(&mut members, left, right)?;

        let right_active = active
            .get_mut(right)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        *right_active = false;

        let left_version = versions
            .get_mut(left)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        *left_version = next_version;

        let left_size_slot = sizes
            .get_mut(left)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        *left_size_slot = merged_size;

        remaining = remaining
            .checked_sub(1)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        for other in 0..size {
            if other == left || active.get(other) != Some(&true) {
                continue;
            }

            let left_distance = pair_distances.get(left, other)?;
            let right_distance = pair_distances.get(right, other)?;

            let updated = updated_linkage_distance(
                linkage,
                left_distance,
                right_distance,
                left_size,
                right_size,
            )?;

            pair_distances.set(left, other, updated)?;

            let (first, second) = ordered_pair(left, other);
            let first_version = versions
                .get(first)
                .copied()
                .ok_or(EnsembleGeometryError::InvalidParameter)?;
            let second_version = versions
                .get(second)
                .copied()
                .ok_or(EnsembleGeometryError::InvalidParameter)?;

            heap.push(MergeCandidate::new(
                updated,
                first,
                second,
                first_version,
                second_version,
            ));
        }
    }

    let clusters = members
        .into_iter()
        .zip(active)
        .filter_map(|(members, active)| active.then_some(members))
        .collect();

    make_clustering(distances, clusters, &[])
}
