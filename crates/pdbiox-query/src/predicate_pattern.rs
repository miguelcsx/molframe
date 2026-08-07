//! Numeric and residue-range pattern parsing for column predicates.

use crate::ast::Column;
use pdbiox_core::contract::{AnalysisPolicy, Namespace};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::{ResidueRef, Structure};

#[derive(Clone, Copy)]
pub(super) enum NumericPattern {
    Value(f64),
    Range(f64, f64),
}

impl NumericPattern {
    pub(super) fn parse(text: &str) -> Result<Self, Diagnostic> {
        if let Some((left, right)) = split_range(text) {
            let left = parse_number(left)?;
            let right = parse_number(right)?;
            return Ok(Self::Range(left.min(right), left.max(right)));
        }
        Ok(Self::Value(parse_number(text)?))
    }

    pub(super) fn matches(self, actual: f64) -> bool {
        match self {
            Self::Value(value) => (actual - value).abs() < f64::EPSILON,
            Self::Range(start, end) => actual >= start && actual <= end,
        }
    }
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

    pub(super) fn matches(&self, actual: &ResidueOrdinal) -> bool {
        match self {
            Self::Value(value) => actual == value,
            Self::Range(start, end) => actual >= start && actual <= end,
        }
    }
}

pub(super) fn residue_ordinal(
    residue: ResidueRef<'_>,
    column: Column,
    policy: &AnalysisPolicy,
) -> Option<ResidueOrdinal> {
    let label = matches!(
        (column, policy.identifiers),
        (Column::LabelResidueId, _) | (Column::ResidueId, Namespace::Label)
    );
    let number = if label {
        residue.label_seq_id()
    } else {
        residue.auth_seq_id().or_else(|| residue.label_seq_id())
    }?;
    Some(ResidueOrdinal {
        number,
        insertion: residue.ins_code().map_or_else(|| "".into(), Into::into),
    })
}

pub(super) fn model_membership(
    structure: &Structure,
    universe: &AtomSelection,
    column: Column,
    patterns: &[NumericPattern],
) -> AtomSelection {
    let matched = if column == Column::ModelIndex {
        (0..structure.model_count())
            .any(|model| patterns.iter().any(|pattern| pattern.matches(model as f64)))
    } else {
        structure
            .data()
            .models()
            .filter_map(pdbiox_core::structure::ModelRef::number)
            .any(|number| {
                patterns
                    .iter()
                    .any(|pattern| pattern.matches(f64::from(number)))
            })
    };
    if matched {
        universe.clone()
    } else {
        AtomSelection::Empty
    }
}

pub(super) fn split_range(text: &str) -> Option<(&str, &str)> {
    if let Some(parts) = text.split_once(':') {
        return Some(parts);
    }
    let lower = text.to_ascii_lowercase();
    if let Some(position) = lower.find("to") {
        return Some((&text[..position], &text[position + 2..]));
    }
    text.char_indices()
        .skip(1)
        .find(|(_, character)| *character == '-')
        .map(|(position, _)| (&text[..position], &text[position + 1..]))
}

fn parse_number(text: &str) -> Result<f64, Diagnostic> {
    text.parse::<f64>()
        .map_err(|_| Diagnostic::new(Code::E4002).with_context("value", text))
}

pub(super) fn parse_residue(text: &str) -> Result<ResidueOrdinal, Diagnostic> {
    let split = text
        .char_indices()
        .skip_while(|(position, character)| *position == 0 && matches!(character, '+' | '-'))
        .find(|(_, character)| !character.is_ascii_digit())
        .map_or(text.len(), |(position, _)| position);
    let number = text[..split]
        .parse::<i32>()
        .map_err(|_| Diagnostic::new(Code::E4002).with_context("value", text))?;
    Ok(ResidueOrdinal {
        number,
        insertion: text[split..].into(),
    })
}
