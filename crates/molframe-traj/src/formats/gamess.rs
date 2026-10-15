//! GAMESS-US optimization and surface-scan output trajectories.

use crate::{FrameValue, Timestep};
use std::collections::BTreeMap;

const BOHR_TO_ANGSTROM: f32 = 0.529_177_2;

/// Coordinate-producing GAMESS run type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GamessRunType {
    /// Geometry optimization (`RUNTYP=OPTIMIZE`).
    Optimize,
    /// Potential-energy surface scan (`RUNTYP=SURFACE`).
    Surface,
}

/// Stable atom identity recovered from GAMESS coordinate records.
#[derive(Clone, Debug, PartialEq)]
pub struct GamessAtom {
    /// GAMESS atom label.
    pub name: Box<str>,
    /// Nuclear charge when the output block carries it.
    pub nuclear_charge: Option<f64>,
}

/// One optimization or surface-scan geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct GamessFrame {
    /// Optimization search ordinal, or surface point ordinal.
    pub step: i64,
    /// Electronic energy in hartree when reported.
    pub energy: Option<f64>,
    /// Surface coordinates 1 and 2, when this is a surface scan.
    pub surface_coordinates: Option<[f64; 2]>,
    /// Cartesian coordinates in ångström.
    pub positions: Vec<[f32; 3]>,
}

impl GamessFrame {
    /// Converts to the common trajectory timestep and retains scalar metadata.
    #[must_use]
    pub fn to_timestep(&self, frame: usize) -> Timestep {
        let mut data = BTreeMap::new();
        data.insert("gamess.step".into(), FrameValue::Integer(self.step));
        if let Some(energy) = self.energy {
            data.insert("gamess.energy_hartree".into(), FrameValue::Float(energy));
        }
        if let Some(coordinates) = self.surface_coordinates {
            data.insert(
                "gamess.surface_coordinates".into(),
                FrameValue::Floats(coordinates.to_vec()),
            );
        }
        Timestep {
            frame,
            positions: self.positions.clone(),
            data,
            ..Timestep::default()
        }
    }
}

/// Fixed atom topology and every geometry in a GAMESS output.
#[derive(Clone, Debug, PartialEq)]
pub struct GamessTrajectory {
    /// Detected coordinate-producing run type.
    pub run_type: GamessRunType,
    /// Atom labels and charges shared by all frames.
    pub atoms: Vec<GamessAtom>,
    /// Frames in output order.
    pub frames: Vec<GamessFrame>,
}

/// Malformed or unsupported GAMESS coordinate output.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum GamessError {
    /// Neither supported run type was declared.
    #[error("GAMESS output is not an OPTIMIZE or SURFACE run")]
    UnsupportedRunType,
    /// The total atom count was absent or invalid.
    #[error("GAMESS total atom count is absent or invalid")]
    MissingAtomCount,
    /// A triggered coordinate block was truncated or non-numeric.
    #[error("invalid GAMESS coordinate block near line {line}")]
    InvalidCoordinate {
        /// One-based source line.
        line: usize,
    },
    /// Atom labels or nuclear charges changed between frames.
    #[error("GAMESS atom identity changed between coordinate blocks")]
    TopologyDrift,
    /// No complete coordinate block was found.
    #[error("GAMESS output contains no complete coordinate frame")]
    NoFrames,
}

/// Parses GAMESS-US/Firefly optimization and surface-scan output coordinates.
///
/// Both ångström and bohr optimization headers are accepted and normalized.
///
/// # Errors
///
/// Returns a run-type, atom-count, coordinate, topology-drift or empty error.
pub fn parse_gamess_output(text: &str) -> Result<GamessTrajectory, GamessError> {
    let lines: Vec<_> = text.lines().collect();
    let run_type = detect_run_type(&lines)?;
    let atom_count = detect_atom_count(&lines)?;
    let (atoms, frames) = match run_type {
        GamessRunType::Optimize => parse_optimization(&lines, atom_count)?,
        GamessRunType::Surface => parse_surface(&lines, atom_count)?,
    };
    if frames.is_empty() {
        return Err(GamessError::NoFrames);
    }
    Ok(GamessTrajectory {
        run_type,
        atoms,
        frames,
    })
}

fn detect_run_type(lines: &[&str]) -> Result<GamessRunType, GamessError> {
    for line in lines {
        if line.contains("RUNTYP=OPTIMIZE") {
            return Ok(GamessRunType::Optimize);
        }
        if line.contains("RUNTYP=SURFACE") {
            return Ok(GamessRunType::Surface);
        }
    }
    Err(GamessError::UnsupportedRunType)
}

fn detect_atom_count(lines: &[&str]) -> Result<usize, GamessError> {
    lines
        .iter()
        .find_map(|line| {
            line.split_once("TOTAL NUMBER OF ATOMS")?
                .1
                .split_once('=')?
                .1
                .split_whitespace()
                .next()?
                .parse()
                .ok()
        })
        .filter(|count| *count > 0)
        .ok_or(GamessError::MissingAtomCount)
}

fn parse_optimization(
    lines: &[&str],
    atom_count: usize,
) -> Result<(Vec<GamessAtom>, Vec<GamessFrame>), GamessError> {
    let energies = optimization_energies(lines);
    let mut topology: Option<Vec<GamessAtom>> = None;
    let mut frames = Vec::new();
    let mut cursor = 0;
    while cursor < lines.len() {
        if lines[cursor].contains("ENERGY=") {
            cursor += 1;
            continue;
        }
        let trigger = lines[cursor]
            .trim_start()
            .trim_start_matches(|character: char| character.is_ascii_digit());
        let Some(step) = trigger
            .trim_start()
            .strip_prefix("NSERCH=")
            .and_then(|tail| tail.split_whitespace().next())
            .and_then(|value| value.parse().ok())
        else {
            cursor += 1;
            continue;
        };
        let Some(header_offset) = lines[cursor + 1..]
            .iter()
            .position(|line| line.contains("COORDINATES OF ALL ATOMS ARE"))
        else {
            cursor += 1;
            continue;
        };
        let header = cursor + 1 + header_offset;
        let scale = if lines[header].contains("BOHR") {
            BOHR_TO_ANGSTROM
        } else {
            1.0
        };
        let Some(separator_offset) = lines[header + 1..]
            .iter()
            .position(|line| line.trim().starts_with('-'))
        else {
            return Err(GamessError::InvalidCoordinate { line: header + 1 });
        };
        let start = header + 2 + separator_offset;
        let (atoms, positions) = read_optimization_atoms(lines, start, atom_count, scale)?;
        accept_topology(&mut topology, atoms)?;
        frames.push(GamessFrame {
            step,
            energy: energies.get(&step).copied(),
            surface_coordinates: None,
            positions,
        });
        cursor = start + atom_count;
    }
    Ok((topology_or_empty(topology), frames))
}

fn optimization_energies(lines: &[&str]) -> BTreeMap<i64, f64> {
    let mut values = BTreeMap::new();
    for line in lines {
        if !line.contains("NSERCH=") || !line.contains("ENERGY=") {
            continue;
        }
        let step = value_after(line, "NSERCH=").and_then(|value| value.parse().ok());
        let energy = value_after(line, "ENERGY=").and_then(|value| value.parse().ok());
        if let (Some(step), Some(energy)) = (step, energy) {
            values.insert(step, energy);
        }
    }
    values
}

fn read_optimization_atoms(
    lines: &[&str],
    start: usize,
    count: usize,
    scale: f32,
) -> Result<(Vec<GamessAtom>, Vec<[f32; 3]>), GamessError> {
    let mut atoms = Vec::with_capacity(count);
    let mut positions = Vec::with_capacity(count);
    for index in 0..count {
        let line_number = start + index + 1;
        let fields: Vec<_> = lines
            .get(start + index)
            .ok_or(GamessError::InvalidCoordinate { line: line_number })?
            .split_whitespace()
            .collect();
        if fields.len() < 5 {
            return Err(GamessError::InvalidCoordinate { line: line_number });
        }
        atoms.push(GamessAtom {
            name: fields[0].into(),
            nuclear_charge: Some(number(fields[1], line_number)?),
        });
        positions.push([
            coordinate(fields[2], scale, line_number)?,
            coordinate(fields[3], scale, line_number)?,
            coordinate(fields[4], scale, line_number)?,
        ]);
    }
    Ok((atoms, positions))
}

fn parse_surface(
    lines: &[&str],
    atom_count: usize,
) -> Result<(Vec<GamessAtom>, Vec<GamessFrame>), GamessError> {
    let mut topology: Option<Vec<GamessAtom>> = None;
    let mut frames = Vec::new();
    let mut cursor = 0;
    while cursor < lines.len() {
        let line = lines[cursor].trim();
        if !line.starts_with("COORD 1=") {
            cursor += 1;
            continue;
        }
        let first = value_after(line, "COORD 1=")
            .and_then(|value| value.parse().ok())
            .ok_or(GamessError::InvalidCoordinate { line: cursor + 1 })?;
        let second = value_after(line, "COORD 2=")
            .and_then(|value| value.parse().ok())
            .ok_or(GamessError::InvalidCoordinate { line: cursor + 1 })?;
        let energy_line = cursor + 1;
        let energy = lines
            .get(energy_line)
            .and_then(|line| line.trim().strip_prefix("HAS ENERGY VALUE"))
            .and_then(|value| value.trim().parse().ok())
            .ok_or(GamessError::InvalidCoordinate {
                line: energy_line + 1,
            })?;
        let start = cursor + 2;
        let (atoms, positions) = read_surface_atoms(lines, start, atom_count)?;
        accept_topology(&mut topology, atoms)?;
        let step = i64::try_from(frames.len())
            .map_err(|_| GamessError::InvalidCoordinate { line: cursor + 1 })?;
        frames.push(GamessFrame {
            step,
            energy: Some(energy),
            surface_coordinates: Some([first, second]),
            positions,
        });
        cursor = start + atom_count;
    }
    Ok((topology_or_empty(topology), frames))
}

fn topology_or_empty(topology: Option<Vec<GamessAtom>>) -> Vec<GamessAtom> {
    let Some(atoms) = topology else {
        return Vec::new();
    };
    atoms
}

fn read_surface_atoms(
    lines: &[&str],
    start: usize,
    count: usize,
) -> Result<(Vec<GamessAtom>, Vec<[f32; 3]>), GamessError> {
    let mut atoms = Vec::with_capacity(count);
    let mut positions = Vec::with_capacity(count);
    for index in 0..count {
        let line_number = start + index + 1;
        let fields: Vec<_> = lines
            .get(start + index)
            .ok_or(GamessError::InvalidCoordinate { line: line_number })?
            .split_whitespace()
            .collect();
        if fields.len() < 4 {
            return Err(GamessError::InvalidCoordinate { line: line_number });
        }
        atoms.push(GamessAtom {
            name: fields[0].into(),
            nuclear_charge: None,
        });
        positions.push([
            coordinate(fields[1], 1.0, line_number)?,
            coordinate(fields[2], 1.0, line_number)?,
            coordinate(fields[3], 1.0, line_number)?,
        ]);
    }
    Ok((atoms, positions))
}

fn accept_topology(
    topology: &mut Option<Vec<GamessAtom>>,
    atoms: Vec<GamessAtom>,
) -> Result<(), GamessError> {
    if topology.as_ref().is_some_and(|expected| expected != &atoms) {
        return Err(GamessError::TopologyDrift);
    }
    if topology.is_none() {
        *topology = Some(atoms);
    }
    Ok(())
}

fn value_after<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    line.split_once(key)?.1.split_whitespace().next()
}

fn number(value: &str, line: usize) -> Result<f64, GamessError> {
    value
        .parse()
        .ok()
        .filter(|value: &f64| value.is_finite())
        .ok_or(GamessError::InvalidCoordinate { line })
}

fn coordinate(value: &str, scale: f32, line: usize) -> Result<f32, GamessError> {
    value
        .parse::<f32>()
        .ok()
        .filter(|value| value.is_finite())
        .map(|value| value * scale)
        .ok_or(GamessError::InvalidCoordinate { line })
}

#[cfg(test)]
#[path = "gamess_tests.rs"]
mod tests;
