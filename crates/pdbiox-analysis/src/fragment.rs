//! Mapping ordered coordinate traces onto explicit fragment libraries.

/// One named fragment reference geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentReference {
    /// Caller-visible stable identifier.
    pub id: Box<str>,
    /// Ordered reference coordinates.
    pub coordinates: Vec<[f32; 3]>,
}

/// Best accepted reference for one complete trace window.
#[derive(Clone, Debug, PartialEq)]
pub struct FragmentMatch {
    /// Window start in the input trace.
    pub start: usize,
    /// Reference identifier.
    pub fragment_id: Box<str>,
    /// Fitted RMSD to the reference.
    pub rmsd: f64,
}

/// Why fragment mapping could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum FragmentMappingError {
    /// Library must be non-empty, with unique IDs and one non-degenerate width.
    #[error("fragment library must be non-empty, uniquely named, finite, and equal-width")]
    InvalidLibrary,
    /// Trace coordinates and acceptance threshold must be finite and threshold non-negative.
    #[error("trace and maximum RMSD must be finite and maximum RMSD non-negative")]
    InvalidInput,
    /// A complete window and reference do not define a unique rigid fit.
    #[error("fragment superposition failed: {0:?}")]
    Superpose(pdbiox_geom::SuperposeError),
}

/// Maps every complete trace window to its closest accepted library fragment.
///
/// Missing trace positions skip only windows that contain them. Equal-RMSD ties
/// choose the lexicographically smaller fragment ID, independent of library order.
///
/// # Errors
///
/// Returns [`FragmentMappingError`] for an invalid library, trace, or threshold.
pub fn map_fragments(
    trace: &[Option<[f32; 3]>],
    library: &[FragmentReference],
    maximum_rmsd: f64,
) -> Result<Vec<FragmentMatch>, FragmentMappingError> {
    let width = validate_library(library)?;
    if !maximum_rmsd.is_finite()
        || maximum_rmsd < 0.0
        || trace
            .iter()
            .flatten()
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(FragmentMappingError::InvalidInput);
    }
    if trace.len() < width {
        return Ok(Vec::new());
    }
    let mut matches = Vec::new();
    for start in 0..=trace.len() - width {
        let Some(window) = complete_window(&trace[start..start + width]) else {
            continue;
        };
        let mut best: Option<(f64, &str)> = None;
        for reference in library {
            let fit = pdbiox_geom::superpose(&window, &reference.coordinates)
                .map_err(FragmentMappingError::Superpose)?;
            let candidate = (fit.rmsd, reference.id.as_ref());
            if best.as_ref().is_none_or(|current| {
                candidate
                    .0
                    .total_cmp(&current.0)
                    .then(candidate.1.cmp(current.1))
                    .is_lt()
            }) {
                best = Some(candidate);
            }
        }
        if let Some((rmsd, id)) = best.filter(|(rmsd, _)| *rmsd <= maximum_rmsd) {
            matches.push(FragmentMatch {
                start,
                fragment_id: id.into(),
                rmsd,
            });
        }
    }
    Ok(matches)
}

fn validate_library(library: &[FragmentReference]) -> Result<usize, FragmentMappingError> {
    let Some(width) = library.first().map(|fragment| fragment.coordinates.len()) else {
        return Err(FragmentMappingError::InvalidLibrary);
    };
    let mut identifiers: Vec<&str> = library
        .iter()
        .map(|fragment| fragment.id.as_ref())
        .collect();
    identifiers.sort_unstable();
    let duplicate = identifiers.windows(2).any(|pair| pair[0] == pair[1]);
    if width < 3
        || duplicate
        || library.iter().any(|fragment| {
            fragment.id.is_empty()
                || fragment.coordinates.len() != width
                || fragment
                    .coordinates
                    .iter()
                    .flatten()
                    .any(|value| !value.is_finite())
        })
    {
        return Err(FragmentMappingError::InvalidLibrary);
    }
    Ok(width)
}

fn complete_window(window: &[Option<[f32; 3]>]) -> Option<Vec<[f32; 3]>> {
    window.iter().copied().collect()
}

#[cfg(test)]
#[path = "fragment_tests.rs"]
mod tests;
