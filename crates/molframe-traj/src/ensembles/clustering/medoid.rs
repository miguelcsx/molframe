// Representative member of a group: the one with least total distance to the rest.
//
// Complexity: O(n^2) validation of the matrix, then O(m^2) time and O(m) space for m members.

/// Returns the member minimizing total within-group distance.
///
/// # Errors
///
/// Returns [`EnsembleGeometryError`] for invalid matrix or empty/out-of-range members.
pub fn medoid(
    distances: &EnsembleDistanceMatrix,
    members: &[usize],
) -> Result<usize, EnsembleGeometryError> {
    validate_matrix(distances)?;
    medoid_validated(distances, members)
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
