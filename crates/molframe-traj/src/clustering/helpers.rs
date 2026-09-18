/// Validates matrix dimensions, finite distances, non-negativity, and symmetry.
fn validate_matrix(distances: &EnsembleDistanceMatrix) -> Result<(), EnsembleGeometryError> {
    let size = distances.size;
    let expected = size
        .checked_mul(size)
        .ok_or(EnsembleGeometryError::DimensionMismatch)?;

    if size == 0 || distances.values.len() != expected {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }

    for row in 0..size {
        let row_values = matrix_row(distances, row)?;
        let diagonal = row_values
            .get(row)
            .copied()
            .ok_or(EnsembleGeometryError::DimensionMismatch)?;

        if !valid_distance(diagonal) {
            return Err(EnsembleGeometryError::InvalidParameter);
        }

        for column in row + 1..size {
            let value = row_values
                .get(column)
                .copied()
                .ok_or(EnsembleGeometryError::DimensionMismatch)?;

            let reverse = matrix_value(distances, column, row)?;

            if !valid_distance(value)
                || !valid_distance(reverse)
                || value.to_bits() != reverse.to_bits()
            {
                return Err(EnsembleGeometryError::InvalidParameter);
            }
        }
    }

    Ok(())
}

/// Returns whether one matrix entry is a valid non-negative distance.
#[inline]
fn valid_distance(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

/// Returns one complete row from the dense distance matrix.
#[inline]
fn matrix_row(
    distances: &EnsembleDistanceMatrix,
    row: usize,
) -> Result<&[f64], EnsembleGeometryError> {
    let start = row
        .checked_mul(distances.size)
        .ok_or(EnsembleGeometryError::DimensionMismatch)?;
    let end = start
        .checked_add(distances.size)
        .ok_or(EnsembleGeometryError::DimensionMismatch)?;

    distances
        .values
        .get(start..end)
        .ok_or(EnsembleGeometryError::DimensionMismatch)
}

/// Returns one checked entry from the dense distance matrix.
#[inline]
fn matrix_value(
    distances: &EnsembleDistanceMatrix,
    row: usize,
    column: usize,
) -> Result<f64, EnsembleGeometryError> {
    let row = matrix_row(distances, row)?;

    row.get(column)
        .copied()
        .ok_or(EnsembleGeometryError::DimensionMismatch)
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

/// Computes every cluster medoid without revalidating the same matrix.
fn cluster_medoids(
    distances: &EnsembleDistanceMatrix,
    members: &[Vec<usize>],
) -> Result<Vec<usize>, EnsembleGeometryError> {
    members
        .iter()
        .map(|cluster| medoid_validated(distances, cluster))
        .collect()
}

/// Computes a medoid after the caller has already validated the matrix.
///
/// Symmetry lets each off-diagonal distance contribute to both candidate totals
/// from a single matrix read while preserving member-order accumulation.
fn medoid_validated(
    distances: &EnsembleDistanceMatrix,
    members: &[usize],
) -> Result<usize, EnsembleGeometryError> {
    if members.is_empty() || members.iter().any(|member| *member >= distances.size) {
        return Err(EnsembleGeometryError::InvalidParameter);
    }

    let mut totals = vec![0.0f64; members.len()];

    for left_position in 0..members.len() {
        let left_member = members[left_position];
        let row = matrix_row(distances, left_member)?;

        let diagonal = row
            .get(left_member)
            .copied()
            .ok_or(EnsembleGeometryError::DimensionMismatch)?;
        totals[left_position] += diagonal;

        for right_position in left_position + 1..members.len() {
            let right_member = members[right_position];
            let distance = row
                .get(right_member)
                .copied()
                .ok_or(EnsembleGeometryError::DimensionMismatch)?;

            totals[left_position] += distance;
            totals[right_position] += distance;
        }
    }

    let mut best: Option<(f64, usize)> = None;

    for (&candidate, &total) in members.iter().zip(&totals) {
        if best.is_none_or(|(best_total, best_candidate)| {
            total
                .total_cmp(&best_total)
                .then(candidate.cmp(&best_candidate))
                == Ordering::Less
        }) {
            best = Some((total, candidate));
        }
    }

    match best {
        Some((_, candidate)) => Ok(candidate),
        None => Err(EnsembleGeometryError::InvalidParameter),
    }
}

/// Collects every observation within `epsilon` of `point` into reusable storage.
fn collect_neighbourhood(
    distances: &EnsembleDistanceMatrix,
    point: usize,
    epsilon: f64,
    output: &mut Vec<usize>,
) -> Result<(), EnsembleGeometryError> {
    output.clear();

    let row = matrix_row(distances, point)?;

    for (other, &distance) in row.iter().enumerate() {
        if distance <= epsilon {
            output.push(other);
        }
    }

    Ok(())
}

/// Creates a reusable index buffer sized for one complete matrix row.
fn reusable_index_buffer(size: usize) -> Result<Vec<usize>, EnsembleGeometryError> {
    let mut buffer = Vec::new();
    buffer
        .try_reserve_exact(size)
        .map_err(|_| EnsembleGeometryError::InvalidParameter)?;
    Ok(buffer)
}

/// State required to expand one DBSCAN cluster.
struct ClusterExpansion<'a> {
    distances: &'a EnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
    label: usize,
    visited: &'a mut [bool],
    labels: &'a mut [Option<usize>],
    queued_epoch: &'a mut [usize],
}

impl ClusterExpansion<'_> {
    /// Expands one core neighbourhood and assigns previously unlabeled points.
    fn run(
        &mut self,
        neighbours: &[usize],
        queue: &mut Vec<usize>,
        scratch: &mut Vec<usize>,
    ) -> Result<(), EnsembleGeometryError> {
        queue.clear();

        for &point in neighbours {
            self.enqueue(point, queue)?;
        }

        let mut head = 0usize;

        while let Some(&point) = queue.get(head) {
            head += 1;

            let was_visited = self
                .visited
                .get(point)
                .copied()
                .ok_or(EnsembleGeometryError::InvalidParameter)?;

            if !was_visited {
                let visited = self
                    .visited
                    .get_mut(point)
                    .ok_or(EnsembleGeometryError::InvalidParameter)?;
                *visited = true;

                collect_neighbourhood(self.distances, point, self.epsilon, scratch)?;

                if scratch.len() >= self.minimum_points {
                    for &neighbour in scratch.iter() {
                        self.enqueue(neighbour, queue)?;
                    }
                }
            }

            let label = self
                .labels
                .get_mut(point)
                .ok_or(EnsembleGeometryError::InvalidParameter)?;

            if label.is_none() {
                *label = Some(self.label);
            }
        }

        Ok(())
    }

    /// Adds a point once to the current cluster's expansion queue.
    fn enqueue(
        &mut self,
        point: usize,
        queue: &mut Vec<usize>,
    ) -> Result<(), EnsembleGeometryError> {
        let epoch = self
            .queued_epoch
            .get_mut(point)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        if *epoch == self.label {
            return Ok(());
        }

        *epoch = self.label;
        queue.push(point);
        Ok(())
    }
}

/// Returns the number of unordered pairs representable for `size` observations.
fn checked_pair_count(size: usize) -> Option<usize> {
    size.checked_mul(size.checked_sub(1)?)?.checked_div(2)
}

/// Orders an unordered pair from the lexicographically smaller identifier.
#[inline]
fn ordered_pair(left: usize, right: usize) -> (usize, usize) {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

/// Normalizes signed zero while retaining every other floating representation.
#[inline]
fn canonical_distance(distance: f64) -> f64 {
    if distance == 0.0 { 0.0 } else { distance }
}
