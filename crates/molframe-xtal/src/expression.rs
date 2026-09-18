//! Assembly operator-expression grammar and bounded Cartesian products.
//!
//! Parsing and expansion are linear in the expression plus the identifiers
//! produced. The product count is checked before any expanded combination is
//! allocated, so hostile ranges cannot turn into unbounded memory use.

use molframe_core::{Code, Diagnostic};
use std::sync::Arc;

/// Default maximum number of generated assembly instances.
pub const DEFAULT_INSTANCE_LIMIT: usize = 100_000;

/// A product of ordered operator-identifier sets.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OperExpression {
    factors: Arc<[Arc<[Box<str>]>]>,
    combinations: usize,
}

impl OperExpression {
    /// Parses an expression under the default assembly-instance limit.
    ///
    /// # Errors
    ///
    /// Returns a registered diagnostic for malformed syntax, invalid ranges or
    /// a Cartesian product larger than the limit.
    pub fn parse(expression: &str) -> Result<Self, Diagnostic> {
        Self::parse_with_limit(expression, DEFAULT_INSTANCE_LIMIT)
    }

    /// Parses an expression under an explicit product limit.
    ///
    /// # Errors
    ///
    /// Returns a registered diagnostic for malformed syntax, invalid ranges or
    /// a Cartesian product larger than `limit`.
    pub fn parse_with_limit(expression: &str, limit: usize) -> Result<Self, Diagnostic> {
        let text = expression.trim();
        if text.is_empty() {
            return Err(malformed(expression));
        }
        let raw_factors = if text.starts_with('(') {
            parenthesised_factors(text).ok_or_else(|| malformed(expression))?
        } else if text.contains(['(', ')']) {
            return Err(malformed(expression));
        } else {
            vec![text]
        };
        let mut factors = Vec::with_capacity(raw_factors.len());
        let mut combinations = 1usize;
        for factor in raw_factors {
            let identifiers = parse_factor(factor).ok_or_else(|| malformed(expression))?;
            combinations = combinations
                .checked_mul(identifiers.len())
                .filter(|count| *count <= limit)
                .ok_or_else(|| {
                    Diagnostic::new(Code::E6011)
                        .with_context("expression", expression)
                        .with_context("limit", limit.to_string())
                })?;
            factors.push(Arc::from(identifiers));
        }
        Ok(Self {
            factors: factors.into(),
            combinations,
        })
    }

    /// Ordered product factors as operator identifiers.
    #[must_use]
    pub fn factors(&self) -> &[Arc<[Box<str>]>] {
        &self.factors
    }

    /// Exact number of combinations the Cartesian product produces.
    #[must_use]
    pub const fn combination_count(&self) -> usize {
        self.combinations
    }

    /// Visits each identifier combination in lexical product order.
    ///
    /// The same fixed-capacity work vector is reused for every callback. The
    /// callback must copy identifiers it needs to retain.
    pub fn for_each_combination(&self, mut visit: impl FnMut(&[&str])) {
        let mut combination = vec![""; self.factors.len()];
        self.visit_factor(0, &mut combination, &mut visit);
    }

    fn visit_factor<'a>(
        &'a self,
        factor: usize,
        combination: &mut [&'a str],
        visit: &mut impl FnMut(&[&str]),
    ) {
        if factor == self.factors.len() {
            visit(combination);
            return;
        }
        for identifier in self.factors[factor].iter() {
            combination[factor] = identifier;
            self.visit_factor(factor + 1, combination, visit);
        }
    }
}

fn parenthesised_factors(text: &str) -> Option<Vec<&str>> {
    let mut factors = Vec::new();
    let bytes = text.as_bytes();
    let mut position = 0usize;
    while position < bytes.len() {
        while bytes.get(position).is_some_and(u8::is_ascii_whitespace) {
            position += 1;
        }
        if bytes.get(position) != Some(&b'(') {
            return None;
        }
        let start = position + 1;
        let relative_end = bytes.get(start..)?.iter().position(|byte| *byte == b')')?;
        let end = start + relative_end;
        if bytes[start..end].contains(&b'(') || text[start..end].trim().is_empty() {
            return None;
        }
        factors.push(text[start..end].trim());
        position = end + 1;
    }
    (!factors.is_empty()).then_some(factors)
}

fn parse_factor(text: &str) -> Option<Vec<Box<str>>> {
    let mut identifiers = Vec::new();
    for part in text.split(',') {
        let item = part.trim();
        if item.is_empty() || item.contains(['(', ')']) {
            return None;
        }
        match numeric_range(item) {
            RangeParse::Valid(start, end) => {
                identifiers.extend((start..=end).map(|value| value.to_string().into_boxed_str()));
            }
            RangeParse::NotRange if valid_identifier(item) => identifiers.push(item.into()),
            RangeParse::NotRange | RangeParse::Invalid => return None,
        }
    }
    (!identifiers.is_empty()).then_some(identifiers)
}

enum RangeParse {
    NotRange,
    Valid(u32, u32),
    Invalid,
}

fn numeric_range(text: &str) -> RangeParse {
    let Some((start, end)) = text.split_once('-') else {
        return RangeParse::NotRange;
    };
    if start.is_empty() || end.is_empty() || end.contains('-') {
        return RangeParse::NotRange;
    }
    let (Ok(start), Ok(end)) = (start.parse::<u32>(), end.parse::<u32>()) else {
        return RangeParse::NotRange;
    };
    if start <= end {
        RangeParse::Valid(start, end)
    } else {
        RangeParse::Invalid
    }
}

fn valid_identifier(text: &str) -> bool {
    text.bytes()
        .all(|byte| byte.is_ascii_graphic() && !matches!(byte, b'(' | b')' | b',' | b'\'' | b'"'))
}

fn malformed(expression: &str) -> Diagnostic {
    Diagnostic::new(Code::E6010).with_context("expression", expression)
}

#[cfg(test)]
#[path = "expression_tests.rs"]
mod tests;
