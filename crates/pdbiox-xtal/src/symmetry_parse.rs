//! Small exact parser helpers for crystallographic coordinate expressions.

use super::{Rational, symmetry_capacity, symmetry_error};
use crate::numeric::usize_to_u32;
use pdbiox_core::Diagnostic;

pub(super) fn parse_component(
    expression: &str,
    coefficients: &mut [i32; 3],
    translation: &mut Rational,
) -> Result<(), Diagnostic> {
    let compact: String = expression
        .chars()
        .filter(|character| !character.is_ascii_whitespace())
        .collect();
    if compact.is_empty() {
        return Err(symmetry_error(expression));
    }
    for (sign, term) in terms(&compact)? {
        if let Some((axis, coefficient)) = variable_term(term)? {
            coefficients[axis] = coefficients[axis]
                .checked_add(sign * coefficient)
                .ok_or_else(symmetry_capacity)?;
        } else {
            let value = parse_rational(term)?;
            *translation =
                translation.checked_add(if sign < 0 { value.negated()? } else { value })?;
        }
    }
    Ok(())
}

fn terms(expression: &str) -> Result<Vec<(i32, &str)>, Diagnostic> {
    let bytes = expression.as_bytes();
    let mut output = Vec::new();
    let mut start = 0usize;
    while start < bytes.len() {
        let (sign, term_start) = match bytes[start] {
            b'+' => (1, start + 1),
            b'-' => (-1, start + 1),
            _ => (1, start),
        };
        let end = bytes[term_start..]
            .iter()
            .position(|byte| matches!(byte, b'+' | b'-'))
            .map_or(bytes.len(), |position| term_start + position);
        let term = expression
            .get(term_start..end)
            .ok_or_else(symmetry_capacity)?;
        if term.is_empty() {
            return Err(symmetry_error(expression));
        }
        output.push((sign, term));
        start = end;
    }
    Ok(output)
}

fn variable_term(term: &str) -> Result<Option<(usize, i32)>, Diagnostic> {
    let Some(variable) = term.chars().last() else {
        return Ok(None);
    };
    let axis = match variable.to_ascii_lowercase() {
        'x' => 0,
        'y' => 1,
        'z' => 2,
        _ => return Ok(None),
    };
    let coefficient = term.strip_suffix(variable).ok_or_else(symmetry_capacity)?;
    let coefficient = if coefficient.is_empty() {
        1
    } else {
        coefficient.parse().map_err(|_| symmetry_error(term))?
    };
    Ok(Some((axis, coefficient)))
}

fn parse_rational(term: &str) -> Result<Rational, Diagnostic> {
    if let Some((numerator, denominator)) = term.split_once('/') {
        if denominator.contains('/') {
            return Err(symmetry_error(term));
        }
        return Rational::new(
            numerator.parse().map_err(|_| symmetry_error(term))?,
            denominator.parse().map_err(|_| symmetry_error(term))?,
        );
    }
    if let Some((whole, fractional)) = term.split_once('.') {
        if fractional.is_empty() || !fractional.bytes().all(|byte| byte.is_ascii_digit()) {
            return Err(symmetry_error(term));
        }
        let denominator = 10_u32
            .checked_pow(usize_to_u32(fractional.len()))
            .ok_or_else(symmetry_capacity)?;
        let whole: i32 = if whole.is_empty() {
            0
        } else {
            whole.parse().map_err(|_| symmetry_error(term))?
        };
        let fraction: u32 = fractional.parse().map_err(|_| symmetry_error(term))?;
        let numerator = i64::from(whole) * i64::from(denominator) + i64::from(fraction);
        return Rational::new(
            i32::try_from(numerator).map_err(|_| symmetry_capacity())?,
            denominator,
        );
    }
    Rational::new(term.parse().map_err(|_| symmetry_error(term))?, 1)
}

pub(super) fn determinant(matrix: [[i32; 3]; 3]) -> i32 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
}

pub(super) fn gcd(mut left: u32, mut right: u32) -> u32 {
    while right != 0 {
        (left, right) = (right, left % right);
    }
    left.max(1)
}
