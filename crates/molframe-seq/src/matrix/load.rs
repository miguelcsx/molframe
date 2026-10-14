//! Resolving a substitution matrix by profile, name or path.
//!
//! Cost: one profile read at most, one `stat` and a hash probe once the table is
//! decoded; building the matrix copies `side^2` scores and a 256-entry lookup.

use super::{ScoreMatrix, SubstitutionMatrix};
use molframe_chem::{TableError, TableRequest, TableStore};

/// Loads the substitution matrix `request` selects: the profile's `score-matrix`
/// line, a cache entry name, a file path, or a table already in hand.
///
/// # Errors
///
/// A [`TableError`] naming how to supply the table when none is selected, or why
/// the selected file is unreadable or malformed.
pub fn load_matrix(
    request: &TableRequest<ScoreMatrix>,
    store: &TableStore,
) -> Result<SubstitutionMatrix, TableError> {
    let provided = request.resolve(store, "matrix")?;
    Ok(SubstitutionMatrix::from_table(&provided))
}

#[cfg(test)]
#[path = "load_tests.rs"]
mod tests;
