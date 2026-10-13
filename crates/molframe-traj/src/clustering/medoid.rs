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
