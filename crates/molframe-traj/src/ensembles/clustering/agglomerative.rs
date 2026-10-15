// Agglomerative clustering to a fixed cluster count over a distance matrix.
//
// Complexity: O(n^2 log n) time and O(n^2) space for n frames, using a lazily
// invalidated merge heap; ties break deterministically.

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

/// Builds the initial heap containing every singleton-cluster pair.
fn initial_merge_heap(
    distances: &PairDistances,
) -> Result<BinaryHeap<MergeCandidate>, EnsembleGeometryError> {
    let pair_count =
        checked_pair_count(distances.size).ok_or(EnsembleGeometryError::InvalidParameter)?;

    let mut heap = BinaryHeap::new();
    heap.try_reserve(pair_count)
        .map_err(|_| EnsembleGeometryError::InvalidParameter)?;

    for left in 0..distances.size {
        for right in left + 1..distances.size {
            heap.push(MergeCandidate::new(
                distances.get(left, right)?,
                left,
                right,
                0,
                0,
            ));
        }
    }

    Ok(heap)
}

/// Removes stale heap entries until a candidate matches current cluster versions.
fn pop_current_candidate(
    heap: &mut BinaryHeap<MergeCandidate>,
    active: &[bool],
    versions: &[usize],
) -> Option<MergeCandidate> {
    while let Some(candidate) = heap.pop() {
        if active.get(candidate.left) != Some(&true) || active.get(candidate.right) != Some(&true) {
            continue;
        }

        if versions.get(candidate.left) != Some(&candidate.left_version)
            || versions.get(candidate.right) != Some(&candidate.right_version)
        {
            continue;
        }

        return Some(candidate);
    }

    None
}

/// Updates one inter-cluster distance after merging two clusters.
fn updated_linkage_distance(
    linkage: Linkage,
    left_distance: f64,
    right_distance: f64,
    left_size: usize,
    right_size: usize,
) -> Result<f64, EnsembleGeometryError> {
    match linkage {
        Linkage::Single => Ok(left_distance.min(right_distance)),
        Linkage::Complete => Ok(left_distance.max(right_distance)),
        Linkage::Average => {
            let left_weight =
                f64_from_usize(left_size).ok_or(EnsembleGeometryError::InvalidParameter)?;
            let right_weight =
                f64_from_usize(right_size).ok_or(EnsembleGeometryError::InvalidParameter)?;
            let total_weight = left_weight + right_weight;

            Ok((left_distance * left_weight + right_distance * right_weight) / total_weight)
        }
    }
}

/// Merges two sorted cluster member lists while keeping the left cluster identity.
fn merge_member_slots(
    members: &mut [Vec<usize>],
    left: usize,
    right: usize,
) -> Result<(), EnsembleGeometryError> {
    if left >= right {
        return Err(EnsembleGeometryError::InvalidParameter);
    }

    let (before_right, from_right) = members.split_at_mut(right);

    let left_members = before_right
        .get_mut(left)
        .ok_or(EnsembleGeometryError::InvalidParameter)?;
    let right_members = from_right
        .first_mut()
        .ok_or(EnsembleGeometryError::InvalidParameter)?;

    let left_members = std::mem::take(left_members);
    let right_members = std::mem::take(right_members);
    let merged = merge_sorted_members(left_members, right_members)?;

    let target = before_right
        .get_mut(left)
        .ok_or(EnsembleGeometryError::InvalidParameter)?;
    *target = merged;

    Ok(())
}

/// Merges two sorted observation lists without an additional sorting pass.
fn merge_sorted_members(
    left: Vec<usize>,
    right: Vec<usize>,
) -> Result<Vec<usize>, EnsembleGeometryError> {
    let total = left
        .len()
        .checked_add(right.len())
        .ok_or(EnsembleGeometryError::InvalidParameter)?;

    let mut merged = Vec::new();
    merged
        .try_reserve_exact(total)
        .map_err(|_| EnsembleGeometryError::InvalidParameter)?;

    let mut left = left.into_iter().peekable();
    let mut right = right.into_iter().peekable();

    loop {
        match (left.peek().copied(), right.peek().copied()) {
            (Some(left_value), Some(right_value)) => {
                if left_value <= right_value {
                    if let Some(value) = left.next() {
                        merged.push(value);
                    }
                } else if let Some(value) = right.next() {
                    merged.push(value);
                }
            }
            (Some(_), None) => {
                merged.extend(left);
                break;
            }
            (None, Some(_)) => {
                merged.extend(right);
                break;
            }
            (None, None) => break,
        }
    }

    Ok(merged)
}

/// Builds labels and representatives from already validated cluster membership.
fn make_clustering(
    distances: &EnsembleDistanceMatrix,
    members: Vec<Vec<usize>>,
    noise: &[usize],
) -> Result<Clustering, EnsembleGeometryError> {
    let mut labels = vec![None; distances.size];

    for (label, cluster) in members.iter().enumerate() {
        for &point in cluster {
            let slot = labels
                .get_mut(point)
                .ok_or(EnsembleGeometryError::InvalidParameter)?;
            *slot = Some(label);
        }
    }

    for &point in noise {
        let slot = labels
            .get_mut(point)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        *slot = None;
    }

    let medoids = cluster_medoids(distances, &members)?;

    Ok(Clustering {
        labels,
        members,
        medoids,
    })
}
