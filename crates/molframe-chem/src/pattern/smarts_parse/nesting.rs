//! Parenthesis-aware splitting for SMARTS atom expressions.

use crate::smarts::SmartsError;

pub(super) fn split_top_level(
    text: &str,
    delimiter: u8,
    offset: usize,
) -> Result<Vec<(usize, &str)>, SmartsError> {
    let mut result = Vec::new();
    let mut start = 0usize;
    let mut depth = 0usize;
    for (position, byte) in text.bytes().enumerate() {
        match byte {
            b'(' => depth += 1,
            b')' if depth == 0 => {
                return Err(SmartsError::new(
                    offset + position,
                    "unmatched closing parenthesis",
                ));
            }
            b')' => depth -= 1,
            value if value == delimiter && depth == 0 => {
                result.push((start, &text[start..position]));
                start = position + 1;
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(SmartsError::new(
            offset + text.len(),
            "unclosed parenthesis",
        ));
    }
    result.push((start, &text[start..]));
    Ok(result)
}
