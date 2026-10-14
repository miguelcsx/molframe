//! Two-pass import of `AlphaFold` PAE JSON into a [`PaeMatrix`].
//!
//! Both published layouts are accepted: the dense `predicted_aligned_error`
//! matrix and the older `residue1`/`residue2`/`distance` triples, each either
//! at the top level or wrapped in a one-element list. The first pass probes the
//! decimal precision and the matrix side, the second fills a buffer sized
//! exactly once, so the transient cost is the final matrix and nothing else.

use super::scan::Scanner;
use crate::{ConfidenceError, PaeMatrix, PaePolicy, PrecisionProbe};
use molframe_core::numeric::f64_to_usize;

#[derive(Default)]
struct Located {
    matrix: Option<usize>,
    first: Option<usize>,
    second: Option<usize>,
    distance: Option<usize>,
    maximum: Option<f64>,
}

/// Reads an `AlphaFold` PAE JSON document into a budgeted dense matrix.
///
/// # Errors
///
/// Rejects unrecognised layouts, ragged or non-square matrices, incomplete or
/// repeated triples, and matrices beyond the policy budget.
pub fn read_pae_json(text: &[u8], policy: PaePolicy) -> Result<PaeMatrix, ConfidenceError> {
    let located = locate(text)?;
    let matrix = if let Some(at) = located.matrix {
        dense(text, at, policy)?
    } else if let (Some(first), Some(second), Some(distance)) =
        (located.first, located.second, located.distance)
    {
        triples(text, [first, second, distance], policy)?
    } else {
        return Err(
            Scanner::new(text).fail("predicted_aligned_error or residue1/residue2/distance")
        );
    };
    Ok(matrix.with_declared_maximum(located.maximum))
}

fn locate(text: &[u8]) -> Result<Located, ConfidenceError> {
    let mut scan = Scanner::new(text);
    if scan.eat(b'[') && scan.peek() == Some(b']') {
        return Err(scan.fail("a PAE object"));
    }
    scan.require(b'{', "an object")?;
    let mut located = Located::default();
    if scan.eat(b'}') {
        return Ok(located);
    }
    loop {
        let key = scan.string()?;
        scan.require(b':', "a colon")?;
        scan.skip_ws();
        let here = scan.offset();
        match key {
            b"predicted_aligned_error" | b"pae" => located.matrix = Some(here),
            b"residue1" => located.first = Some(here),
            b"residue2" => located.second = Some(here),
            b"distance" => located.distance = Some(here),
            b"max_predicted_aligned_error" => {
                located.maximum = Some(scan.number()?);
                if !scan.eat(b',') {
                    return Ok(located);
                }
                continue;
            }
            _ => {}
        }
        scan.skip_value()?;
        if !scan.eat(b',') {
            return Ok(located);
        }
    }
}

/// Visits every cell of a matrix literal and returns `(rows, columns)`.
fn walk_matrix(
    text: &[u8],
    at: usize,
    mut visit: impl FnMut(usize, usize, f64) -> Result<(), ConfidenceError>,
) -> Result<(usize, usize), ConfidenceError> {
    let mut scan = Scanner::at(text, at);
    scan.require(b'[', "a matrix")?;
    let (mut rows, mut columns) = (0_usize, 0_usize);
    if scan.eat(b']') {
        return Ok((0, 0));
    }
    loop {
        scan.require(b'[', "a matrix row")?;
        let mut length = 0_usize;
        if !scan.eat(b']') {
            loop {
                visit(rows, length, scan.number()?)?;
                length += 1;
                if !scan.eat(b',') {
                    break;
                }
            }
            scan.require(b']', "the end of a row")?;
        }
        if rows == 0 {
            columns = length;
        } else if length != columns {
            return Err(ConfidenceError::Ragged {
                row: rows,
                expected: columns,
                found: length,
            });
        }
        rows += 1;
        if !scan.eat(b',') {
            break;
        }
    }
    scan.require(b']', "the end of the matrix")?;
    Ok((rows, columns))
}

fn dense(text: &[u8], at: usize, policy: PaePolicy) -> Result<PaeMatrix, ConfidenceError> {
    let mut probe = PrecisionProbe::new();
    let (rows, columns) = walk_matrix(text, at, |_, _, value| {
        probe.observe_f64(value);
        Ok(())
    })?;
    if rows != columns {
        return Err(ConfidenceError::SizeMismatch {
            residues: rows,
            matrix: columns,
        });
    }
    let layout = policy.layout(&probe, rows)?;
    let mut matrix = PaeMatrix::zeroed(rows, layout, policy.max_bytes)?;
    walk_matrix(text, at, |row, column, value| {
        matrix.set_f64(row, column, value)
    })?;
    Ok(matrix)
}

/// Visits aligned `(residue1, residue2, distance)` elements of three arrays.
fn walk_triples(
    text: &[u8],
    at: [usize; 3],
    mut visit: impl FnMut(usize, usize, f64) -> Result<(), ConfidenceError>,
) -> Result<(), ConfidenceError> {
    let mut scans = at.map(|position| Scanner::at(text, position));
    for scan in &mut scans {
        scan.require(b'[', "an array")?;
    }
    if scans[0].eat(b']') {
        let closed = scans[1].eat(b']') && scans[2].eat(b']');
        return if closed {
            Ok(())
        } else {
            Err(scans[1].fail("an empty array"))
        };
    }
    loop {
        let first = residue(&mut scans[0])?;
        let second = residue(&mut scans[1])?;
        let value = scans[2].number()?;
        visit(first, second, value)?;
        let more = [scans[0].eat(b','), scans[1].eat(b','), scans[2].eat(b',')];
        if more == [false; 3] {
            break;
        }
        if more != [true; 3] {
            return Err(scans[0].fail("arrays of equal length"));
        }
    }
    for scan in &mut scans {
        scan.require(b']', "the end of an array")?;
    }
    Ok(())
}

fn residue(scan: &mut Scanner<'_>) -> Result<usize, ConfidenceError> {
    let at = scan.offset();
    let value = scan.number()?;
    if value >= 1.0 && value.fract() == 0.0 {
        Ok(f64_to_usize(value))
    } else {
        Err(Scanner::at(&[], at).fail("a one-based residue number"))
    }
}

fn triples(text: &[u8], at: [usize; 3], policy: PaePolicy) -> Result<PaeMatrix, ConfidenceError> {
    let mut probe = PrecisionProbe::new();
    let mut size = 0_usize;
    walk_triples(text, at, |first, second, value| {
        probe.observe_f64(value);
        size = size.max(first).max(second);
        Ok(())
    })?;
    let layout = policy.layout(&probe, size)?;
    let mut matrix = PaeMatrix::zeroed(size, layout, policy.max_bytes)?;
    let mut seen = vec![0_u64; (size * size).div_ceil(64)];
    let mut count = 0_u64;
    walk_triples(text, at, |first, second, value| {
        let (row, column) = (first - 1, second - 1);
        let cell = row * size + column;
        let (word, bit) = (cell / 64, 1_u64 << (cell % 64));
        if seen[word] & bit != 0 {
            return Err(ConfidenceError::Conflict { row, column });
        }
        seen[word] |= bit;
        count += 1;
        matrix.set_f64(row, column, value)
    })?;
    let expected = (size as u64) * (size as u64);
    if count != expected {
        return Err(ConfidenceError::Incomplete {
            expected,
            found: count,
        });
    }
    Ok(matrix)
}

#[cfg(test)]
#[path = "import_tests.rs"]
mod tests;
