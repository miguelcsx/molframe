//! `OpenDX` scalar grid reader.
//!
//! This is the format electrostatics solvers write a potential map in. Only the
//! regular-lattice subset is accepted — `gridpositions` with three `delta`
//! lines and a rank-zero array — because that is the shape a scalar field on a
//! box actually takes, and quietly mis-reading a general DX object would be
//! worse than refusing it.
//!
//! Geometry is already in angstroms. Values run with the third axis fastest and
//! are transposed once on the way in.

use super::{Axes, GridError, canonical_value_index, density_map, value_count, zeroed_values};
use crate::mrc::DensityMap;

const FORMAT: &str = "dx";

/// Reads an `OpenDX` regular scalar grid.
///
/// # Errors
///
/// Returns [`GridError`] for non-UTF-8 text, a header that is not a regular
/// lattice, dimensions out of range, or a truncated or non-numeric value block.
pub fn read_dx(bytes: &[u8]) -> Result<DensityMap, GridError> {
    let Ok(text) = std::str::from_utf8(bytes) else {
        return Err(GridError::NotUtf8 { format: FORMAT });
    };
    let mut lines = text.lines();
    let header = read_header(&mut lines)?;
    let total = value_count(header.axes.counts, FORMAT)?;
    if header.declared_items.is_some_and(|items| items != total) {
        return Err(GridError::InvalidHeader {
            format: FORMAT,
            reason: "declared item count disagrees with the grid counts",
        });
    }
    let values = read_values(&mut lines, header.axes.counts, total)?;
    Ok(density_map(&header.axes, values))
}

/// Grid geometry and the declared item count, if the file stated one.
struct DxHeader {
    axes: Axes,
    declared_items: Option<usize>,
}

/// Consumes header statements up to and including `data follows`.
fn read_header<'a>(lines: &mut impl Iterator<Item = &'a str>) -> Result<DxHeader, GridError> {
    let mut counts = None;
    let mut origin = None;
    let mut steps = [[0.0_f64; 3]; 3];
    let mut delta_index = 0;
    let mut declared_items = None;

    for line in lines {
        let trimmed = line.trim_start();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("object") {
            if let Some(after) = rest.find("counts") {
                let parsed = triple_usize(&rest[after + "counts".len()..])?;
                // `gridconnections` repeats the same counts; the first
                // statement to carry them is the authority.
                if counts.is_none() {
                    counts = Some(parsed);
                }
            } else if let Some(after) = rest.find("items") {
                declared_items = rest[after + "items".len()..]
                    .split_ascii_whitespace()
                    .next()
                    .and_then(|token| token.parse::<usize>().ok());
                if rest.contains("data follows") {
                    break;
                }
            }
        } else if let Some(rest) = trimmed.strip_prefix("origin") {
            origin = Some(triple_f64(rest)?);
        } else if let Some(rest) = trimmed.strip_prefix("delta") {
            if delta_index >= 3 {
                return Err(GridError::InvalidHeader {
                    format: FORMAT,
                    reason: "more than three delta vectors",
                });
            }
            steps[delta_index] = triple_f64(rest)?;
            delta_index += 1;
        }
    }

    let (Some(counts), Some(origin)) = (counts, origin) else {
        return Err(GridError::InvalidHeader {
            format: FORMAT,
            reason: "a regular grid needs counts and an origin",
        });
    };
    if delta_index != 3 {
        return Err(GridError::InvalidHeader {
            format: FORMAT,
            reason: "a regular grid needs three delta vectors",
        });
    }
    Ok(DxHeader {
        axes: Axes {
            origin,
            steps,
            counts,
        },
        declared_items,
    })
}

/// Reads exactly `total` values, stopping at the statement that closes the
/// array rather than running off into trailing attributes.
fn read_values<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    counts: [usize; 3],
    total: usize,
) -> Result<Vec<f32>, GridError> {
    let mut values = zeroed_values(total, FORMAT)?;
    let mut parsed = 0_usize;
    for line in lines {
        let trimmed = line.trim_start();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed.starts_with("attribute") || trimmed.starts_with("object") {
            break;
        }
        for token in trimmed.split_ascii_whitespace() {
            if parsed == total {
                break;
            }
            match token.parse::<f32>() {
                Ok(value) => {
                    values[canonical_value_index(parsed, counts, 1)] = value;
                    parsed += 1;
                }
                Err(_) => {
                    return Err(GridError::InvalidNumber {
                        format: FORMAT,
                        index: parsed,
                    });
                }
            }
        }
        if parsed == total {
            break;
        }
    }
    if parsed == total {
        Ok(values)
    } else {
        Err(GridError::Truncated { format: FORMAT })
    }
}

fn triple_usize(rest: &str) -> Result<[usize; 3], GridError> {
    let mut fields = rest.split_ascii_whitespace();
    let mut counts = [0_usize; 3];
    for count in &mut counts {
        let Some(token) = fields.next() else {
            return Err(GridError::InvalidHeader {
                format: FORMAT,
                reason: "grid counts",
            });
        };
        match token.parse::<usize>() {
            Ok(value) => *count = value,
            Err(_) => {
                return Err(GridError::InvalidHeader {
                    format: FORMAT,
                    reason: "grid counts",
                });
            }
        }
    }
    Ok(counts)
}

fn triple_f64(rest: &str) -> Result<[f64; 3], GridError> {
    let mut fields = rest.split_ascii_whitespace();
    let mut values = [0.0_f64; 3];
    for value in &mut values {
        let Some(token) = fields.next() else {
            return Err(GridError::InvalidHeader {
                format: FORMAT,
                reason: "a three-component vector",
            });
        };
        match token.parse::<f64>() {
            Ok(parsed) => *value = parsed,
            Err(_) => {
                return Err(GridError::InvalidHeader {
                    format: FORMAT,
                    reason: "a three-component vector",
                });
            }
        }
    }
    Ok(values)
}
