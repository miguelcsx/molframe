//! Bounded similar-word enumeration under any substitution scorer.

use crate::Score;

/// Explicit score and resource limits for similar-kmer enumeration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SimilarKmerOptions {
    /// Minimum sum of substitution scores.
    pub minimum_score: i32,
    /// Maximum results retained after deterministic ranking.
    pub max_results: usize,
    /// Maximum Cartesian candidates the caller permits the search to visit.
    pub max_candidates: usize,
}

/// One candidate word and its summed substitution score.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SimilarKmer {
    /// Candidate word.
    pub word: Vec<u8>,
    /// Sum of position-wise substitution scores.
    pub score: i32,
}

/// Similar-word policy cannot be executed within its declared bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum SimilarKmerError {
    /// Alphabet, word, or result limit is empty.
    EmptyInput,
    /// Alphabet contains duplicate symbols.
    DuplicateSymbol,
    /// Cartesian candidate count exceeds the explicit ceiling.
    CandidateLimit {
        /// Complete Cartesian candidate count.
        required: usize,
        /// Caller-declared visitation ceiling.
        allowed: usize,
    },
    /// A complete candidate score exceeded the exact supported range.
    NumericOverflow,
}

impl std::fmt::Display for SimilarKmerError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => {
                formatter.write_str("similar k-mer input and limits must be nonempty")
            }
            Self::DuplicateSymbol => {
                formatter.write_str("similar k-mer alphabet contains duplicates")
            }
            Self::CandidateLimit { required, allowed } => write!(
                formatter,
                "similar k-mer enumeration needs {required} candidates, limit is {allowed}"
            ),
            Self::NumericOverflow => {
                formatter.write_str("similar k-mer score exceeds the supported exact range")
            }
        }
    }
}

impl std::error::Error for SimilarKmerError {}

/// Enumerates and ranks every bounded word meeting a substitution-score threshold.
///
/// # Errors
///
/// Returns [`SimilarKmerError`] for empty/duplicate inputs, numeric overflow,
/// or when the complete Cartesian search exceeds `max_candidates`; it never
/// publishes a truncated search as complete.
pub fn similar_kmers(
    query: &[u8],
    alphabet: &[u8],
    scorer: &impl Score,
    options: SimilarKmerOptions,
) -> Result<Vec<SimilarKmer>, SimilarKmerError> {
    if query.is_empty() || alphabet.is_empty() || options.max_results == 0 {
        return Err(SimilarKmerError::EmptyInput);
    }
    let mut unique = alphabet.to_vec();
    unique.sort_unstable();
    unique.dedup();
    if unique.len() != alphabet.len() {
        return Err(SimilarKmerError::DuplicateSymbol);
    }
    let required = (0..query.len()).try_fold(1usize, |count, _| count.checked_mul(alphabet.len()));
    let Some(required) = required else {
        return Err(SimilarKmerError::CandidateLimit {
            required: usize::MAX,
            allowed: options.max_candidates,
        });
    };
    if required > options.max_candidates {
        return Err(SimilarKmerError::CandidateLimit {
            required,
            allowed: options.max_candidates,
        });
    }
    let mut search = SimilarSearch {
        query,
        alphabet,
        scorer,
        minimum: i64::from(options.minimum_score),
        word: Vec::with_capacity(query.len()),
        results: Vec::new(),
    };
    search.enumerate(0, 0)?;
    let mut results = search.results;
    results.sort_by(|left, right| {
        right
            .score
            .cmp(&left.score)
            .then(left.word.cmp(&right.word))
    });
    results.truncate(options.max_results);
    Ok(results)
}

struct SimilarSearch<'a, S> {
    query: &'a [u8],
    alphabet: &'a [u8],
    scorer: &'a S,
    minimum: i64,
    word: Vec<u8>,
    results: Vec<SimilarKmer>,
}

impl<S: Score> SimilarSearch<'_, S> {
    fn enumerate(&mut self, position: usize, score: i64) -> Result<(), SimilarKmerError> {
        if position == self.query.len() {
            let score = i32::try_from(score).map_err(|_| SimilarKmerError::NumericOverflow)?;
            if i64::from(score) >= self.minimum {
                self.results.push(SimilarKmer {
                    word: self.word.clone(),
                    score,
                });
            }
            return Ok(());
        }
        for symbol in self.alphabet {
            self.word.push(*symbol);
            let next_position = position
                .checked_add(1)
                .ok_or(SimilarKmerError::NumericOverflow)?;
            let next_score = score
                .checked_add(i64::from(self.scorer.score(self.query[position], *symbol)))
                .ok_or(SimilarKmerError::NumericOverflow)?;
            self.enumerate(next_position, next_score)?;
            self.word.pop();
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "similar_tests.rs"]
mod tests;
