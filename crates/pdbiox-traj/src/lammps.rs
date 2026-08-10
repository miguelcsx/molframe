//! LAMMPS native text dump reader with custom-column discovery.

use crate::numeric::f32_triplet;
use crate::{FrameValue, Timestep};
use std::collections::{BTreeMap, BTreeSet};

/// Malformed or unsupported LAMMPS dump content.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum LammpsError {
    /// A required `ITEM:` record was absent or out of order.
    #[error("expected LAMMPS record {expected}")]
    MissingRecord {
        /// Required header prefix.
        expected: &'static str,
    },
    /// A numeric field could not be parsed.
    #[error("invalid LAMMPS numeric field {field}")]
    InvalidNumber {
        /// Column or record name.
        field: Box<str>,
    },
    /// No supported Cartesian or scaled coordinate triplet was declared.
    #[error("LAMMPS atom columns do not contain x/y/z, xu/yu/zu or xs/ys/zs")]
    MissingCoordinates,
    /// An atom row has fewer values than its header.
    #[error("LAMMPS atom row is shorter than the declared column list")]
    ShortAtomRow,
    /// Atom identifiers repeat within a frame.
    #[error("duplicate LAMMPS atom id {id}")]
    DuplicateAtomId {
        /// Repeated identifier.
        id: i64,
    },
    /// Frames do not carry a stable atom count.
    #[error("LAMMPS frame atom count changed from {expected} to {found}")]
    AtomCountMismatch {
        /// Count established by the first frame.
        expected: usize,
        /// Count in the current frame.
        found: usize,
    },
}

/// Parses native `dump atom` and `dump custom` snapshots.
///
/// Atom columns are discovered from `ITEM: ATOMS`. Rows are sorted by `id` when
/// present, preventing processor-dependent dump order from changing topology.
/// Wrapped, un-wrapped and scaled coordinate triplets are supported; velocity
/// and force triplets populate the typed auxiliary buffers.
///
/// # Errors
///
/// Returns a structural error for missing records, invalid numbers, duplicate
/// atom IDs, changing atom counts or absent coordinate columns.
pub fn parse_lammps_dump(text: &str) -> Result<Vec<Timestep>, LammpsError> {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let mut frames = Vec::new();
    let mut expected_atoms = None;
    while let Some(header) = lines.next() {
        require(header, "ITEM: TIMESTEP")?;
        let step = integer(lines.next(), "timestep")?;
        require(
            next(&mut lines, "ITEM: NUMBER OF ATOMS")?,
            "ITEM: NUMBER OF ATOMS",
        )?;
        let atom_count = unsigned(lines.next(), "atom count")?;
        if let Some(expected) = expected_atoms
            && expected != atom_count
        {
            return Err(LammpsError::AtomCountMismatch {
                expected,
                found: atom_count,
            });
        }
        expected_atoms = Some(atom_count);
        let bounds_header = next(&mut lines, "ITEM: BOX BOUNDS")?;
        require(bounds_header, "ITEM: BOX BOUNDS")?;
        let triclinic = bounds_header.split_whitespace().any(|field| field == "xy");
        let bounds = read_bounds(&mut lines, triclinic)?;
        let atom_header = next(&mut lines, "ITEM: ATOMS")?;
        require(atom_header, "ITEM: ATOMS")?;
        let columns: Vec<_> = atom_header.split_whitespace().skip(2).collect();
        let layout = Layout::new(&columns)?;
        let mut rows = Vec::with_capacity(atom_count);
        for ordinal in 0..atom_count {
            let fields: Vec<_> = next(&mut lines, "atom row")?.split_whitespace().collect();
            if fields.len() < columns.len() {
                return Err(LammpsError::ShortAtomRow);
            }
            rows.push(layout.read(&fields, ordinal, &bounds)?);
        }
        if layout.id.is_some() {
            rows.sort_unstable_by_key(|row| row.id);
            let mut ids = BTreeSet::new();
            for row in &rows {
                if !ids.insert(row.id) {
                    return Err(LammpsError::DuplicateAtomId { id: row.id });
                }
            }
        }
        frames.push(to_timestep(frames.len(), step, &rows, &bounds));
    }
    Ok(frames)
}

fn require(found: &str, expected: &'static str) -> Result<(), LammpsError> {
    if found.starts_with(expected) {
        Ok(())
    } else {
        Err(LammpsError::MissingRecord { expected })
    }
}

fn next<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    expected: &'static str,
) -> Result<&'a str, LammpsError> {
    lines.next().ok_or(LammpsError::MissingRecord { expected })
}

fn integer(line: Option<&str>, field: &'static str) -> Result<i64, LammpsError> {
    line.and_then(|value| value.trim().parse().ok())
        .ok_or_else(|| LammpsError::InvalidNumber {
            field: field.into(),
        })
}

fn unsigned(line: Option<&str>, field: &'static str) -> Result<usize, LammpsError> {
    line.and_then(|value| value.trim().parse().ok())
        .ok_or_else(|| LammpsError::InvalidNumber {
            field: field.into(),
        })
}

#[derive(Clone, Copy)]
struct Bounds {
    origin: [f64; 3],
    basis: [[f64; 3]; 3],
}

fn read_bounds<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    triclinic: bool,
) -> Result<Bounds, LammpsError> {
    let mut raw = [[0.0; 3]; 3];
    for (axis, row) in raw.iter_mut().enumerate() {
        let values: Vec<_> = next(lines, "box bound")?
            .split_whitespace()
            .map(str::parse::<f64>)
            .collect();
        if values.len() < 2 || values.iter().any(Result::is_err) {
            return Err(LammpsError::InvalidNumber {
                field: format!("box bound {axis}").into(),
            });
        }
        for (slot, value) in row.iter_mut().zip(values.into_iter()) {
            *slot = value.map_err(|_| LammpsError::InvalidNumber {
                field: format!("box bound {axis}").into(),
            })?;
        }
    }
    let [xy, xz, yz] = if triclinic {
        [raw[0][2], raw[1][2], raw[2][2]]
    } else {
        [0.0; 3]
    };
    let xlo = raw[0][0] - 0.0_f64.min(xy).min(xz).min(xy + xz);
    let xhi = raw[0][1] - 0.0_f64.max(xy).max(xz).max(xy + xz);
    let ylo = raw[1][0] - 0.0_f64.min(yz);
    let yhi = raw[1][1] - 0.0_f64.max(yz);
    let zlo = raw[2][0];
    let zhi = raw[2][1];
    Ok(Bounds {
        origin: [xlo, ylo, zlo],
        basis: [
            [xhi - xlo, xy, xz],
            [0.0, yhi - ylo, yz],
            [0.0, 0.0, zhi - zlo],
        ],
    })
}

#[derive(Clone, Copy)]
enum CoordinateKind {
    Cartesian,
    Scaled,
}

struct Layout {
    id: Option<usize>,
    coordinates: [usize; 3],
    kind: CoordinateKind,
    velocity: Option<[usize; 3]>,
    force: Option<[usize; 3]>,
}

impl Layout {
    fn new(columns: &[&str]) -> Result<Self, LammpsError> {
        let find = |name| columns.iter().position(|column| *column == name);
        let (coordinates, kind) = triplet(columns, ["xu", "yu", "zu"])
            .map(|indices| (indices, CoordinateKind::Cartesian))
            .or_else(|| {
                triplet(columns, ["x", "y", "z"])
                    .map(|indices| (indices, CoordinateKind::Cartesian))
            })
            .or_else(|| {
                triplet(columns, ["xs", "ys", "zs"])
                    .map(|indices| (indices, CoordinateKind::Scaled))
            })
            .ok_or(LammpsError::MissingCoordinates)?;
        Ok(Self {
            id: find("id"),
            coordinates,
            kind,
            velocity: triplet(columns, ["vx", "vy", "vz"]),
            force: triplet(columns, ["fx", "fy", "fz"]),
        })
    }

    fn read(
        &self,
        fields: &[&str],
        ordinal: usize,
        bounds: &Bounds,
    ) -> Result<AtomRow, LammpsError> {
        let ordinal = i64::try_from(ordinal).map_err(|_| LammpsError::InvalidNumber {
            field: "atom ordinal".into(),
        })?;
        let id = self
            .id
            .map_or(Ok(ordinal), |index| parse(fields[index], "id"))?;
        let coordinates = values(fields, self.coordinates, "coordinates")?;
        let position = match self.kind {
            CoordinateKind::Cartesian => narrow(coordinates, "coordinates")?,
            CoordinateKind::Scaled => scaled(coordinates, bounds)?,
        };
        Ok(AtomRow {
            id,
            position,
            velocity: self
                .velocity
                .map(|indices| values(fields, indices, "velocity"))
                .transpose()?
                .map(|values| narrow(values, "velocity"))
                .transpose()?,
            force: self
                .force
                .map(|indices| values(fields, indices, "force"))
                .transpose()?
                .map(|values| narrow(values, "force"))
                .transpose()?,
        })
    }
}

fn triplet(columns: &[&str], names: [&str; 3]) -> Option<[usize; 3]> {
    let find = |name| columns.iter().position(|column| *column == name);
    Some([find(names[0])?, find(names[1])?, find(names[2])?])
}

fn parse(value: &str, field: &'static str) -> Result<i64, LammpsError> {
    value.parse().map_err(|_| LammpsError::InvalidNumber {
        field: field.into(),
    })
}

fn values(
    fields: &[&str],
    indices: [usize; 3],
    field: &'static str,
) -> Result<[f64; 3], LammpsError> {
    let mut output = [0.0; 3];
    for (slot, index) in output.iter_mut().zip(indices) {
        *slot = fields[index]
            .parse()
            .map_err(|_| LammpsError::InvalidNumber {
                field: field.into(),
            })?;
    }
    Ok(output)
}

fn scaled(fractional: [f64; 3], bounds: &Bounds) -> Result<[f32; 3], LammpsError> {
    narrow(
        [
            bounds.origin[0]
                + bounds.basis[0][0] * fractional[0]
                + bounds.basis[0][1] * fractional[1]
                + bounds.basis[0][2] * fractional[2],
            bounds.origin[1]
                + bounds.basis[1][1] * fractional[1]
                + bounds.basis[1][2] * fractional[2],
            bounds.origin[2] + bounds.basis[2][2] * fractional[2],
        ],
        "scaled coordinates",
    )
}

fn narrow(values: [f64; 3], field: &'static str) -> Result<[f32; 3], LammpsError> {
    f32_triplet(values).ok_or_else(|| LammpsError::InvalidNumber {
        field: field.into(),
    })
}

struct AtomRow {
    id: i64,
    position: [f32; 3],
    velocity: Option<[f32; 3]>,
    force: Option<[f32; 3]>,
}

fn to_timestep(index: usize, step: i64, rows: &[AtomRow], bounds: &Bounds) -> Timestep {
    let velocities = rows.iter().map(|row| row.velocity).collect();
    let forces = rows.iter().map(|row| row.force).collect();
    let mut data = BTreeMap::new();
    data.insert("lammps.step".into(), FrameValue::Integer(step));
    data.insert(
        "lammps.origin".into(),
        FrameValue::Floats(bounds.origin.to_vec()),
    );
    Timestep {
        frame: index,
        positions: rows.iter().map(|row| row.position).collect(),
        velocities,
        forces,
        cell: crate::cell::cell_from_vectors([
            [bounds.basis[0][0], 0.0, 0.0],
            [bounds.basis[0][1], bounds.basis[1][1], 0.0],
            [bounds.basis[0][2], bounds.basis[1][2], bounds.basis[2][2]],
        ]),
        data,
        ..Timestep::default()
    }
}

#[cfg(test)]
#[path = "lammps_tests.rs"]
mod tests;
