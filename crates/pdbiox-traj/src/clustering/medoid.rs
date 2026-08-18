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
