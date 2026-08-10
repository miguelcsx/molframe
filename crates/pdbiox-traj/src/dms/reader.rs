//! Strict, single-pass reads of DMS chemical-structure tables.

use crate::dms::schema::{
    BOND, GLOBAL_CELL, PARTICLE, PARTICLE_FLOATS, PARTICLE_INTEGERS, PARTICLE_PROJECTION,
    PARTICLE_TEXT, REQUIRED_BOND, REQUIRED_CELL, REQUIRED_PARTICLE, VERSION,
    optional_table_columns, require_columns, select_projection, table_columns,
    validate_standard_types,
};
use crate::dms::{
    DmsBond, DmsCell, DmsError, DmsFrame, DmsParticle, DmsSystem, DmsTopology, DmsVersion,
};
use rusqlite::{Connection, OpenFlags, Row};
use std::path::Path;

#[cfg(test)]
use std::collections::BTreeSet;

/// Reads topology, coordinates, velocities, and the periodic cell from a DMS database.
///
/// Particle identifiers and bond references are validated before the model is
/// returned. Optional standard particle columns retain SQL `NULL` as `None`.
///
/// # Errors
///
/// Returns an error when `SQLite` cannot read the file, a structural table is
/// malformed, identifiers violate DMS ordering, or numeric values are invalid.
pub fn read_dms(path: impl AsRef<Path>) -> Result<DmsSystem, DmsError> {
    let connection = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )?;
    let version = read_version(&connection)?;
    let (particles, positions, velocities) = read_particles(&connection)?;
    let bonds = read_bonds(&connection, particles.len())?;
    let cell = read_cell(&connection)?;
    Ok(DmsSystem {
        version,
        topology: DmsTopology { particles, bonds },
        frame: DmsFrame {
            positions,
            velocities,
            cell,
        },
    })
}

type ParticleRows = (Vec<DmsParticle>, Vec<[f64; 3]>, Vec<Option<[f64; 3]>>);
type ParticleRecord = (DmsParticle, [f64; 3], Option<[f64; 3]>);

fn read_particles(connection: &Connection) -> Result<ParticleRows, DmsError> {
    let columns = table_columns(connection, PARTICLE)?;
    require_columns(PARTICLE, &columns, &REQUIRED_PARTICLE)?;
    validate_standard_types(
        connection,
        PARTICLE,
        &PARTICLE_INTEGERS,
        &PARTICLE_FLOATS,
        &PARTICLE_TEXT,
    )?;
    let sql = format!(
        "SELECT {} FROM {PARTICLE} ORDER BY id",
        select_projection(&columns, &PARTICLE_PROJECTION)
    );
    let mut statement = connection.prepare(&sql)?;
    let mut rows = statement.query([])?;
    let mut particles = Vec::new();
    let mut positions = Vec::new();
    let mut velocities = Vec::new();
    while let Some(row) = rows.next()? {
        let id: i64 = row.get(0)?;
        let expected = i64::try_from(particles.len()).map_err(|_| DmsError::InvalidParticleIds)?;
        if id != expected {
            return Err(DmsError::InvalidParticleIds);
        }
        let (particle, position, velocity) = particle(row, id)?;
        particles.push(particle);
        positions.push(position);
        velocities.push(velocity);
    }
    Ok((particles, positions, velocities))
}

fn particle(row: &Row<'_>, id: i64) -> Result<ParticleRecord, DmsError> {
    let position = [row.get(4)?, row.get(5)?, row.get(6)?];
    if position.iter().any(|value: &f64| !value.is_finite()) {
        return Err(DmsError::InvalidParticle { id });
    }
    let velocity = vector(row.get(7)?, row.get(8)?, row.get(9)?, id)?;
    let atomic_number = optional_u16(row.get(1)?, id)?;
    let particle = DmsParticle {
        atomic_number,
        component: row.get(2)?,
        nonbonded_type: row.get(3)?,
        mass: finite(row.get(10)?, id)?,
        charge: finite(row.get(11)?, id)?,
        residue_id: row.get(12)?,
        residue_name: text(row.get(13)?),
        chain: text(row.get(14)?),
        segment: text(row.get(15)?),
        name: text(row.get(16)?),
        insertion: text(row.get(17)?),
        formal_charge: finite(row.get(18)?, id)?,
        occupancy: finite(row.get(19)?, id)?,
        b_factor: finite(row.get(20)?, id)?,
        temperature_group: row.get(21)?,
        energy_group: row.get(22)?,
        ligand_group: row.get(23)?,
        bias_group: row.get(24)?,
    };
    Ok((particle, position, velocity))
}

fn vector(
    x: Option<f64>,
    y: Option<f64>,
    z: Option<f64>,
    id: i64,
) -> Result<Option<[f64; 3]>, DmsError> {
    match (x, y, z) {
        (None, None, None) => Ok(None),
        (Some(x), Some(y), Some(z)) if [x, y, z].iter().all(|value| value.is_finite()) => {
            Ok(Some([x, y, z]))
        }
        _ => Err(DmsError::InvalidParticle { id }),
    }
}

fn optional_u16(value: Option<i64>, id: i64) -> Result<Option<u16>, DmsError> {
    value
        .map(|number| u16::try_from(number).map_err(|_| DmsError::InvalidParticle { id }))
        .transpose()
}

fn finite(value: Option<f64>, id: i64) -> Result<Option<f64>, DmsError> {
    match value {
        Some(number) if number.is_finite() => Ok(Some(number)),
        Some(_) => Err(DmsError::InvalidParticle { id }),
        None => Ok(None),
    }
}

fn text(value: Option<String>) -> Option<Box<str>> {
    value.map(String::into_boxed_str)
}

fn read_bonds(connection: &Connection, particles: usize) -> Result<Vec<DmsBond>, DmsError> {
    let columns = table_columns(connection, BOND)?;
    require_columns(BOND, &columns, &REQUIRED_BOND)?;
    validate_standard_types(connection, BOND, &["p0", "p1"], &["order"], &[])?;
    let mut statement = connection.prepare("SELECT p0, p1, \"order\" FROM bond ORDER BY p0, p1")?;
    let rows = statement.query_map([], |row| {
        Ok((
            row.get::<_, i64>(0)?,
            row.get::<_, i64>(1)?,
            row.get::<_, f64>(2)?,
        ))
    })?;
    let mut bonds = Vec::new();
    let mut previous = None;
    for row in rows {
        let (p0, p1, order) = row?;
        let pair = (p0, p1);
        if p0 < 0
            || p0 >= p1
            || usize::try_from(p1).map_or(true, |index| index >= particles)
            || !order.is_finite()
            || previous == Some(pair)
        {
            return Err(DmsError::InvalidBond { p0, p1 });
        }
        bonds.push(DmsBond {
            p0: u32::try_from(p0).map_err(|_| DmsError::InvalidBond { p0, p1 })?,
            p1: u32::try_from(p1).map_err(|_| DmsError::InvalidBond { p0, p1 })?,
            order,
        });
        previous = Some(pair);
    }
    Ok(bonds)
}

fn read_cell(connection: &Connection) -> Result<Option<DmsCell>, DmsError> {
    let Some(columns) = optional_table_columns(connection, GLOBAL_CELL)? else {
        return Ok(None);
    };
    require_columns(GLOBAL_CELL, &columns, &REQUIRED_CELL)?;
    validate_standard_types(connection, GLOBAL_CELL, &["id"], &["x", "y", "z"], &[])?;
    let mut statement = connection.prepare("SELECT id, x, y, z FROM global_cell ORDER BY id")?;
    let rows = statement.query_map([], |row| {
        Ok((
            row.get::<_, i64>(0)?,
            [row.get(1)?, row.get(2)?, row.get(3)?],
        ))
    })?;
    let collected = rows.collect::<Result<Vec<(i64, [f64; 3])>, _>>()?;
    if collected.len() != 3 {
        return Err(DmsError::InvalidCell);
    }
    let mut vectors = [[0.0; 3]; 3];
    for (expected, (id, vector)) in collected.into_iter().enumerate() {
        let expected_id = i64::try_from(expected).map_err(|_| DmsError::InvalidCell)?;
        if id != expected_id || vector.iter().any(|value| !value.is_finite()) {
            return Err(DmsError::InvalidCell);
        }
        vectors[expected] = vector;
    }
    Ok(Some(DmsCell { vectors }))
}

fn read_version(connection: &Connection) -> Result<Option<DmsVersion>, DmsError> {
    let Some(columns) = optional_table_columns(connection, VERSION)? else {
        return Ok(None);
    };
    require_columns(VERSION, &columns, &["major", "minor"])?;
    validate_standard_types(connection, VERSION, &["major", "minor"], &[], &[])?;
    let mut statement = connection.prepare("SELECT major, minor FROM dms_version")?;
    let rows = statement
        .query_map([], |row| Ok((row.get::<_, i64>(0)?, row.get::<_, i64>(1)?)))?
        .collect::<Result<Vec<_>, _>>()?;
    let [(major, minor)] = rows.as_slice() else {
        return Err(DmsError::InvalidVersion);
    };
    Ok(Some(DmsVersion {
        major: u32::try_from(*major).map_err(|_| DmsError::InvalidVersion)?,
        minor: u32::try_from(*minor).map_err(|_| DmsError::InvalidVersion)?,
    }))
}

#[cfg(test)]
pub(crate) fn column_set(names: &[&str]) -> BTreeSet<Box<str>> {
    names.iter().map(|name| Box::<str>::from(*name)).collect()
}
