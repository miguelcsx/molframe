//! GROMOS11 block trajectory (`TRC`/`TRJ`) reader.

use crate::numeric::f32_triplet;
use crate::{FrameValue, Timestep};
use molframe_core::structure::UnitCell;
use std::collections::BTreeMap;

const NM_TO_ANGSTROM: f32 = 10.0;

/// Boundary convention declared by a GROMOS `GENBOX` block.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GromosBoundary {
    /// No periodic boundary conditions (`0`).
    Vacuum,
    /// Rectangular periodic cell (`1`).
    Rectangular,
    /// Triclinic periodic cell (`2`).
    Triclinic,
    /// Truncated-octahedral convention (`-1`).
    TruncatedOctahedron,
}

/// Complete block trajectory and per-frame boundary semantics.
#[derive(Clone, Debug, PartialEq)]
pub struct GromosTrajectory {
    /// Mandatory `TITLE` block content.
    pub title: Box<str>,
    /// Frames in source order.
    pub frames: Vec<Timestep>,
    /// Boundary convention corresponding to each frame.
    pub boundaries: Vec<Option<GromosBoundary>>,
}

/// Invalid GROMOS11 trajectory text.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum GromosError {
    /// A named block is missing its terminating `END` line.
    #[error("unterminated GROMOS11 block: {name}")]
    Unterminated {
        /// Block name.
        name: Box<str>,
    },
    /// Mandatory title or coordinate blocks are absent.
    #[error("GROMOS11 trajectory is missing {block}")]
    Missing {
        /// Required block name.
        block: &'static str,
    },
    /// A block has malformed or non-finite values.
    #[error("invalid GROMOS11 {block} block")]
    Invalid {
        /// Invalid block name.
        block: &'static str,
    },
    /// Coordinate counts changed between frames.
    #[error("GROMOS11 atom count changed from {expected} to {actual}")]
    AtomCount {
        /// First-frame atom count.
        expected: usize,
        /// Current-frame atom count.
        actual: usize,
    },
}

/// Parses standard `TIMESTEP`, `POSITIONRED` and `GENBOX` blocks.
///
/// Unknown blocks are retained as textual `Timestep::data` entries under the
/// `gromos11.<BLOCK>` key instead of being silently discarded. Shifted origins,
/// Euler rotations and truncated-octahedral boundary identity are retained too.
///
/// # Errors
///
/// Returns the first unterminated, missing, malformed or topology-drift error.
pub fn parse_gromos11_trc(source: &str) -> Result<GromosTrajectory, GromosError> {
    let blocks = blocks(source)?;
    let title = blocks
        .iter()
        .find(|block| block.name.as_ref() == "TITLE")
        .map(|block| block.lines.join("\n").trim().into())
        .ok_or(GromosError::Missing { block: "TITLE" })?;
    let mut drafts = Vec::new();
    let mut current = Draft::default();
    for block in blocks
        .into_iter()
        .filter(|block| block.name.as_ref() != "TITLE")
    {
        match block.name.as_ref() {
            "TIMESTEP" => {
                if current.positions.is_some() {
                    drafts.push(current);
                    current = Draft::default();
                }
                current.timestep = Some(parse_timestep(&block.lines)?);
            }
            "POSITIONRED" => {
                if current.positions.is_some() {
                    drafts.push(current);
                    current = Draft::default();
                }
                current.positions = Some(parse_positions(&block.lines)?);
            }
            "GENBOX" => current.cell = Some(parse_genbox(&block.lines)?),
            name => {
                current.unknown.insert(
                    format!("gromos11.{name}").into(),
                    FrameValue::Text(block.lines.join("\n").into()),
                );
            }
        }
    }
    if current.positions.is_some() {
        drafts.push(current);
    }
    if drafts.is_empty() {
        return Err(GromosError::Missing {
            block: "POSITIONRED",
        });
    }
    build_trajectory(title, drafts)
}

#[derive(Default)]
struct Draft {
    timestep: Option<(i64, f64)>,
    positions: Option<Vec<[f32; 3]>>,
    cell: Option<BoxData>,
    unknown: BTreeMap<Box<str>, FrameValue>,
}

struct BoxData {
    boundary: GromosBoundary,
    cell: Option<UnitCell>,
    origin: [f64; 3],
    euler: [f64; 3],
}

fn build_trajectory(title: Box<str>, drafts: Vec<Draft>) -> Result<GromosTrajectory, GromosError> {
    let expected = drafts[0].positions.as_ref().map_or(0, Vec::len);
    let mut frames = Vec::with_capacity(drafts.len());
    let mut boundaries = Vec::with_capacity(drafts.len());
    for (ordinal, draft) in drafts.into_iter().enumerate() {
        let positions = draft.positions.ok_or(GromosError::Missing {
            block: "POSITIONRED",
        })?;
        if positions.len() != expected {
            return Err(GromosError::AtomCount {
                expected,
                actual: positions.len(),
            });
        }
        let mut data = draft.unknown;
        let time = draft.timestep.map(|(_, time)| time);
        if let Some((step, _)) = draft.timestep {
            data.insert("trajectory_step".into(), FrameValue::Integer(step));
        }
        let (boundary, cell) = match draft.cell {
            Some(box_data) => {
                data.insert(
                    "gromos11.origin_nm".into(),
                    FrameValue::Floats(box_data.origin.into()),
                );
                data.insert(
                    "gromos11.euler_degrees".into(),
                    FrameValue::Floats(box_data.euler.into()),
                );
                (Some(box_data.boundary), box_data.cell)
            }
            None => (None, None),
        };
        frames.push(Timestep {
            frame: ordinal,
            time,
            positions,
            cell,
            data,
            ..Timestep::default()
        });
        boundaries.push(boundary);
    }
    populate_dt(&mut frames);
    Ok(GromosTrajectory {
        title,
        frames,
        boundaries,
    })
}

struct Block {
    name: Box<str>,
    lines: Vec<Box<str>>,
}

fn blocks(source: &str) -> Result<Vec<Block>, GromosError> {
    let lines: Vec<&str> = source.lines().collect();
    let mut cursor = 0;
    let mut result = Vec::new();
    while cursor < lines.len() {
        let name = lines[cursor].trim();
        cursor += 1;
        if name.is_empty() || name.starts_with('#') {
            continue;
        }
        let mut content = Vec::new();
        let mut terminated = false;
        while cursor < lines.len() {
            let line = lines[cursor];
            cursor += 1;
            if line.trim() == "END" {
                terminated = true;
                break;
            }
            content.push(line.into());
        }
        if !terminated {
            return Err(GromosError::Unterminated { name: name.into() });
        }
        result.push(Block {
            name: name.into(),
            lines: content,
        });
    }
    Ok(result)
}

fn parse_timestep(lines: &[Box<str>]) -> Result<(i64, f64), GromosError> {
    let values: Vec<&str> = data_lines(lines).flat_map(str::split_whitespace).collect();
    if values.len() != 2 {
        return Err(invalid("TIMESTEP"));
    }
    let step = values[0].parse().map_err(|_| invalid("TIMESTEP"))?;
    let time = values[1].parse().map_err(|_| invalid("TIMESTEP"))?;
    if !f64::is_finite(time) {
        return Err(invalid("TIMESTEP"));
    }
    Ok((step, time))
}

fn parse_positions(lines: &[Box<str>]) -> Result<Vec<[f32; 3]>, GromosError> {
    let mut positions = Vec::new();
    for line in data_lines(lines) {
        let values = parse_three(line, "POSITIONRED")?;
        positions.push(
            f32_triplet(values.map(|value| value * f64::from(NM_TO_ANGSTROM)))
                .ok_or_else(|| invalid("POSITIONRED"))?,
        );
    }
    if positions.is_empty()
        || positions
            .iter()
            .flatten()
            .any(|coordinate| !coordinate.is_finite())
    {
        return Err(invalid("POSITIONRED"));
    }
    Ok(positions)
}

fn parse_genbox(lines: &[Box<str>]) -> Result<BoxData, GromosError> {
    let content: Vec<&str> = data_lines(lines).collect();
    let boundary = match content
        .first()
        .and_then(|line| line.trim().parse::<i32>().ok())
    {
        Some(0) => GromosBoundary::Vacuum,
        Some(1) => GromosBoundary::Rectangular,
        Some(2) => GromosBoundary::Triclinic,
        Some(-1) => GromosBoundary::TruncatedOctahedron,
        _ => return Err(invalid("GENBOX")),
    };
    if boundary == GromosBoundary::Vacuum {
        return Ok(BoxData {
            boundary,
            cell: None,
            origin: [0.0; 3],
            euler: [0.0; 3],
        });
    }
    if content.len() != 5 {
        return Err(invalid("GENBOX"));
    }
    let lengths = parse_three(content[1], "GENBOX")?.map(|value| value * 10.0);
    let angles = parse_three(content[2], "GENBOX")?;
    let origin = parse_three(content[3], "GENBOX")?;
    let euler = parse_three(content[4], "GENBOX")?;
    let cell = UnitCell { lengths, angles };
    if lengths
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        || angles.iter().any(|value| !value.is_finite())
        || origin.iter().chain(&euler).any(|value| !value.is_finite())
    {
        return Err(invalid("GENBOX"));
    }
    Ok(BoxData {
        boundary,
        cell: Some(cell),
        origin,
        euler,
    })
}

fn data_lines(lines: &[Box<str>]) -> impl Iterator<Item = &str> {
    lines
        .iter()
        .map(AsRef::as_ref)
        .filter(|line: &&str| !line.trim().is_empty() && !line.trim_start().starts_with('#'))
}

fn parse_three(line: &str, block: &'static str) -> Result<[f64; 3], GromosError> {
    let mut values = line.split_whitespace().map(str::parse::<f64>);
    let result = [
        values
            .next()
            .ok_or_else(|| invalid(block))?
            .map_err(|_| invalid(block))?,
        values
            .next()
            .ok_or_else(|| invalid(block))?
            .map_err(|_| invalid(block))?,
        values
            .next()
            .ok_or_else(|| invalid(block))?
            .map_err(|_| invalid(block))?,
    ];
    if values.next().is_some() {
        Err(invalid(block))
    } else {
        Ok(result)
    }
}

fn populate_dt(frames: &mut [Timestep]) {
    for index in 1..frames.len() {
        if let (Some(previous), Some(current)) = (frames[index - 1].time, frames[index].time) {
            frames[index].dt = Some(current - previous);
        }
    }
    if frames.len() > 1 {
        frames[0].dt = frames[1].dt;
    }
}

const fn invalid(block: &'static str) -> GromosError {
    GromosError::Invalid { block }
}

#[cfg(test)]
#[path = "trc_tests.rs"]
mod tests;
