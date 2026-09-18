//! Numeric and residue-range pattern parsing for column predicates.

use crate::ast::Column;
use molframe_core::contract::{AnalysisPolicy, Namespace};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::ResidueRef;
use std::cmp::Ordering;

#[derive(Clone, Copy)]
pub(super) enum NumericPattern {
    Value(f64),
    Range(f64, f64),
}

impl NumericPattern {
    /// Parses a scalar or ascending-normalized inclusive numeric range.
    pub(super) fn parse(text: &str) -> Result<Self, Diagnostic> {
        if let Some((left, right)) = split_range(text) {
            let left = parse_number(left)?;
            let right = parse_number(right)?;
            return Ok(Self::Range(left.min(right), left.max(right)));
        }
        Ok(Self::Value(parse_number(text)?))
    }

    /// Returns in `O(1)` whether `actual` satisfies this value or range.
    #[inline]
    pub(super) fn matches(self, actual: f64) -> bool {
        match self {
            Self::Value(value) => (actual - value).abs() < f64::EPSILON,
            Self::Range(start, end) => actual >= start && actual <= end,
        }
    }
}

/// Compiled numeric patterns supporting logarithmic membership checks.
pub(super) struct NumericMatcher {
    kind: NumericMatcherKind,
}

/// Internal representation whose invariants can only be built by `NumericMatcher`.
enum NumericMatcherKind {
    Empty,
    Single(NumericPattern),
    Compiled {
        values: Vec<f64>,
        ranges: Vec<(f64, f64)>,
    },
}

impl NumericMatcher {
    /// Compiles numeric patterns for repeated matching.
    ///
    /// Zero- and one-pattern inputs allocate no memory. Larger inputs are
    /// normalized in `O(P log P)` time and `O(P)` space, where `P` is the
    /// number of patterns.
    pub(super) fn from_patterns(patterns: &[NumericPattern]) -> Self {
        let kind = match patterns {
            [] => NumericMatcherKind::Empty,
            [pattern] => NumericMatcherKind::Single(*pattern),
            _ => {
                let value_count = patterns
                    .iter()
                    .filter(|pattern| matches!(pattern, NumericPattern::Value(_)))
                    .count();
                let mut values = Vec::with_capacity(value_count);
                let range_count = match patterns.len().checked_sub(value_count) {
                    Some(count) => count,
                    None => 0,
                };
                let mut ranges = Vec::with_capacity(range_count);

                for &pattern in patterns {
                    match pattern {
                        NumericPattern::Value(value) if value.is_finite() => {
                            values.push(value);
                        }
                        NumericPattern::Range(start, end)
                            if !start.is_nan() && !end.is_nan() && start <= end =>
                        {
                            ranges.push((start, end));
                        }
                        _ => {}
                    }
                }

                values.sort_unstable_by(f64::total_cmp);
                values.dedup_by(|left, right| left.total_cmp(right) == Ordering::Equal);
                merge_numeric_ranges(&mut ranges);

                NumericMatcherKind::Compiled { values, ranges }
            }
        };

        Self { kind }
    }

    /// Returns whether `actual` matches any compiled pattern.
    ///
    /// Matching is `O(1)` for zero or one pattern and `O(log P)` for a
    /// compiled multi-pattern matcher, with no per-call allocation.
    #[inline]
    pub(super) fn matches(&self, actual: f64) -> bool {
        match &self.kind {
            NumericMatcherKind::Empty => false,
            NumericMatcherKind::Single(pattern) => pattern.matches(actual),
            NumericMatcherKind::Compiled { values, ranges } => {
                matches_numeric_range(ranges, actual) || matches_numeric_value(values, actual)
            }
        }
    }
}

/// Merges overlapping inclusive numeric ranges in place.
///
/// The input is sorted as part of the operation. Runtime is `O(R log R)` and
/// additional space is `O(1)` beyond the vector itself.
fn merge_numeric_ranges(ranges: &mut Vec<(f64, f64)>) {
    if ranges.len() < 2 {
        return;
    }

    ranges.sort_unstable_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
    });

    let mut write = 0usize;
    for read in 1..ranges.len() {
        let (start, end) = ranges[read];
        let current_end = ranges[write].1;

        if start <= current_end {
            if end > current_end {
                ranges[write].1 = end;
            }
        } else {
            write += 1;
            ranges[write] = (start, end);
        }
    }

    ranges.truncate(write + 1);
}

/// Tests a value against sorted, non-overlapping inclusive ranges.
///
/// Runtime is `O(log R)` and no memory is allocated.
#[inline]
fn matches_numeric_range(ranges: &[(f64, f64)], actual: f64) -> bool {
    let position = ranges.partition_point(|(start, _)| *start <= actual);

    position
        .checked_sub(1)
        .and_then(|index| ranges.get(index))
        .is_some_and(|(_, end)| actual <= *end)
}

/// Tests a value against sorted scalar patterns using the original epsilon rule.
///
/// Only the two numeric neighbors around `actual` can satisfy the fixed
/// `f64::EPSILON` threshold, so matching runs in `O(log V)` time without
/// allocation.
#[inline]
fn matches_numeric_value(values: &[f64], actual: f64) -> bool {
    if !actual.is_finite() {
        return false;
    }

    let position = values.partition_point(|value| *value < actual);

    values
        .get(position)
        .is_some_and(|value| (actual - *value).abs() < f64::EPSILON)
        || position
            .checked_sub(1)
            .and_then(|index| values.get(index))
            .is_some_and(|value| (actual - *value).abs() < f64::EPSILON)
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct ResidueOrdinal {
    pub(super) number: i32,
    pub(super) insertion: Box<str>,
}

pub(super) enum ResiduePattern {
    Value(ResidueOrdinal),
    Range(ResidueOrdinal, ResidueOrdinal),
}

impl ResiduePattern {
    /// Parses a residue ordinal or inclusive residue range from `text`.
    ///
    /// Endpoints are normalized according to `ResidueOrdinal` ordering.
    /// Parsing runs in `O(text.len())` time with one owned insertion string per
    /// parsed endpoint.
    pub(super) fn parse(text: &str) -> Result<Self, Diagnostic> {
        if let Some((left, right)) = split_range(text) {
            let left = parse_residue(left)?;
            let right = parse_residue(right)?;

            return Ok(if left <= right {
                Self::Range(left, right)
            } else {
                Self::Range(right, left)
            });
        }

        Ok(Self::Value(parse_residue(text)?))
    }

    /// Matches borrowed residue components without constructing a
    /// `ResidueOrdinal`.
    ///
    /// Runtime is `O(1)` and no heap allocation is performed.
    #[inline]
    pub(super) fn matches_parts(&self, number: i32, insertion: &str) -> bool {
        match self {
            Self::Value(value) => number == value.number && insertion == value.insertion.as_ref(),
            Self::Range(start, end) => {
                compare_residue_parts(number, insertion, start) != Ordering::Less
                    && compare_residue_parts(number, insertion, end) != Ordering::Greater
            }
        }
    }

    /// Returns the inclusive lower endpoint of this residue pattern.
    #[inline]
    fn start(&self) -> &ResidueOrdinal {
        match self {
            Self::Value(value) | Self::Range(value, _) => value,
        }
    }

    /// Returns the inclusive upper endpoint of this residue pattern.
    #[inline]
    fn end(&self) -> &ResidueOrdinal {
        match self {
            Self::Value(value) | Self::Range(_, value) => value,
        }
    }
}

/// Compiled residue patterns represented as a single pattern or sorted intervals.
pub(super) struct ResidueMatcher {
    single: Option<ResiduePattern>,
    intervals: Vec<ResidueInterval>,
}

impl ResidueMatcher {
    /// Compiles residue patterns for repeated matching.
    ///
    /// Zero- and one-pattern inputs allocate no additional memory. Larger
    /// inputs are normalized in `O(P log P)` time and `O(P)` space.
    pub(super) fn from_patterns(mut patterns: Vec<ResiduePattern>) -> Self {
        patterns.retain(|pattern| match pattern {
            ResiduePattern::Value(_) => true,
            ResiduePattern::Range(start, end) => start <= end,
        });

        match patterns.len() {
            0 => Self {
                single: None,
                intervals: Vec::new(),
            },
            1 => {
                let mut iterator = patterns.into_iter();

                Self {
                    single: iterator.next(),
                    intervals: Vec::new(),
                }
            }
            _ => {
                patterns.sort_unstable_by(|left, right| {
                    left.start()
                        .cmp(right.start())
                        .then_with(|| left.end().cmp(right.end()))
                });

                let mut intervals = Vec::<ResidueInterval>::with_capacity(patterns.len());

                for pattern in patterns {
                    let interval = ResidueInterval::from_pattern(pattern);

                    if let Some(last) = intervals.last_mut()
                        && interval.start <= *last.end()
                    {
                        let extends = interval.end() > last.end();

                        if extends {
                            last.end = Some(interval.into_end());
                        }

                        continue;
                    }

                    intervals.push(interval);
                }

                Self {
                    single: None,
                    intervals,
                }
            }
        }
    }

    /// Returns whether borrowed residue components match any compiled pattern.
    ///
    /// Matching is `O(1)` for zero or one pattern and `O(log P)` for compiled
    /// interval sets, with no per-call allocation.
    #[inline]
    pub(super) fn matches_parts(&self, number: i32, insertion: &str) -> bool {
        if let Some(pattern) = &self.single {
            return pattern.matches_parts(number, insertion);
        }

        let position = self.intervals.partition_point(|interval| {
            compare_residue_parts(number, insertion, &interval.start) != Ordering::Less
        });

        position
            .checked_sub(1)
            .and_then(|index| self.intervals.get(index))
            .is_some_and(|interval| {
                compare_residue_parts(number, insertion, interval.end()) != Ordering::Greater
            })
    }
}

/// Owned interval used by `ResidueMatcher` after normalization.
struct ResidueInterval {
    start: ResidueOrdinal,
    end: Option<ResidueOrdinal>,
}

impl ResidueInterval {
    /// Converts a parsed residue pattern into an owned interval without cloning.
    fn from_pattern(pattern: ResiduePattern) -> Self {
        match pattern {
            ResiduePattern::Value(start) => Self { start, end: None },
            ResiduePattern::Range(start, end) => Self {
                start,
                end: Some(end),
            },
        }
    }

    /// Returns the inclusive upper endpoint, reusing `start` for point intervals.
    #[inline]
    fn end(&self) -> &ResidueOrdinal {
        match &self.end {
            Some(end) => end,
            None => &self.start,
        }
    }

    /// Consumes the interval and returns its inclusive upper endpoint without cloning.
    fn into_end(self) -> ResidueOrdinal {
        match self.end {
            Some(end) => end,
            None => self.start,
        }
    }
}

/// Compares borrowed residue components against an owned residue ordinal.
///
/// Comparison is lexicographic by residue number and insertion code and
/// allocates no memory.
#[inline]
fn compare_residue_parts(number: i32, insertion: &str, ordinal: &ResidueOrdinal) -> Ordering {
    number
        .cmp(&ordinal.number)
        .then_with(|| insertion.cmp(ordinal.insertion.as_ref()))
}

/// Returns the numeric residue identifier selected by `column` and `policy`.
///
/// The lookup is `O(1)` and performs no allocation.
#[inline]
pub(super) fn residue_number(
    residue: ResidueRef<'_>,
    column: Column,
    policy: &AnalysisPolicy,
) -> Option<i32> {
    let label = matches!(
        (column, policy.identifiers),
        (Column::LabelResidueId, _) | (Column::ResidueId, Namespace::Label)
    );

    if label {
        residue.label_seq_id()
    } else {
        residue.auth_seq_id()
    }
}

/// Returns the residue insertion code as a borrowed string, using `""` when absent.
///
/// The lookup is `O(1)` and performs no allocation.
#[inline]
pub(super) fn residue_insertion(residue: ResidueRef<'_>) -> &str {
    match residue.ins_code() {
        Some(insertion) => insertion,
        None => "",
    }
}

/// Splits scalar/range syntax using `:`, case-insensitive `to`, or `-` delimiters.
///
/// Delimiters retain their original precedence. Runtime is `O(text.len())` and
/// no temporary lowercase string is allocated.
pub(super) fn split_range(text: &str) -> Option<(&str, &str)> {
    if let Some(parts) = text.split_once(':') {
        return Some(parts);
    }

    if let Some(position) = find_ascii_to(text) {
        return Some((&text[..position], &text[position + 2..]));
    }

    text.char_indices()
        .skip(1)
        .find(|(_, character)| *character == '-')
        .map(|(position, _)| (&text[..position], &text[position + 1..]))
}

/// Finds the first ASCII case-insensitive `"to"` delimiter without allocation.
///
/// Runtime is `O(text.len())` and the returned offset is always a UTF-8 boundary.
#[inline]
fn find_ascii_to(text: &str) -> Option<usize> {
    text.as_bytes().windows(2).position(|window| {
        window[0].eq_ignore_ascii_case(&b't') && window[1].eq_ignore_ascii_case(&b'o')
    })
}

/// Parses a floating-point scalar and maps parse failures to `E4002`.
///
/// Runtime is `O(text.len())` with no allocation on the success path.
fn parse_number(text: &str) -> Result<f64, Diagnostic> {
    text.parse::<f64>()
        .map_err(|_| Diagnostic::new(Code::E4002).with_context("value", text))
}

/// Parses a residue number followed by an optional insertion code.
///
/// Runtime is `O(text.len())`; exactly one owned insertion string is created in
/// the returned ordinal.
pub(super) fn parse_residue(text: &str) -> Result<ResidueOrdinal, Diagnostic> {
    let bytes = text.as_bytes();

    let digit_start = match bytes.first() {
        Some(&b'+' | &b'-') => 1,
        _ => 0,
    };

    let split = match bytes[digit_start..]
        .iter()
        .position(|byte| !byte.is_ascii_digit())
    {
        Some(offset) => digit_start + offset,
        None => text.len(),
    };

    let number = text[..split]
        .parse::<i32>()
        .map_err(|_| Diagnostic::new(Code::E4002).with_context("value", text))?;

    Ok(ResidueOrdinal {
        number,
        insertion: text[split..].into(),
    })
}
