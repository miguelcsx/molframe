//! Deterministic DMS writer for the standard chemical-structure schema.

use crate::dms::{DmsBond, DmsError, DmsParticle, DmsSystem};
use rusqlite::{Connection, OpenFlags, Transaction, params};
use std::collections::BTreeSet;
use std::fs::OpenOptions;
use std::io::ErrorKind;
use std::path::Path;

const PARTICLE_SCHEMA: &str = "
CREATE TABLE particle (
 id INTEGER PRIMARY KEY,
 anum INTEGER,
 msys_ct INTEGER,
 nbtype INTEGER,
 x REAL NOT NULL,
 y REAL NOT NULL,
 z REAL NOT NULL,
 vx REAL,
 vy REAL,
 vz REAL,
 mass REAL,
 charge REAL,
 resid INTEGER,
 resname TEXT,
 chain TEXT,
 segid TEXT,
 name TEXT,
 insertion TEXT,
 formal_charge REAL,
 occupancy REAL,
 bfactor REAL,
 grp_temperature INTEGER,
 grp_energy INTEGER,
 grp_ligand INTEGER,
 grp_bias INTEGER
);
CREATE TABLE bond (
 p0 INTEGER NOT NULL,
 p1 INTEGER NOT NULL,
 \"order\" REAL NOT NULL,
 PRIMARY KEY (p0, p1),
 CHECK (p0 < p1),
 FOREIGN KEY (p0) REFERENCES particle(id),
 FOREIGN KEY (p1) REFERENCES particle(id)
);";

const PARTICLE_INSERT: &str = "
INSERT INTO particle (
 id, anum, msys_ct, nbtype, x, y, z, vx, vy, vz, mass, charge,
 resid, resname, chain, segid, name, insertion, formal_charge,
 occupancy, bfactor, grp_temperature, grp_energy, grp_ligand, grp_bias
) VALUES (
 ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12,
 ?13, ?14, ?15, ?16, ?17, ?18, ?19, ?20, ?21, ?22, ?23, ?24, ?25
)";

/// Writes a new DMS `SQLite` database without replacing an existing path.
///
/// Rows are emitted in particle-id and bond-id order. Standard DMS units are
/// retained directly: Angstrom, picosecond, elementary charge, and atomic mass.
///
/// # Errors
///
/// Returns an error for an existing destination, inconsistent arrays, invalid
/// numeric data or connectivity, or any `SQLite` failure.
pub fn write_dms(path: impl AsRef<Path>, system: &DmsSystem) -> Result<(), DmsError> {
    validate(system)?;
    let path = path.as_ref();
    reserve_destination(path)?;
    let mut connection = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_WRITE | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )?;
    connection.execute_batch("PRAGMA foreign_keys = ON")?;
    let transaction = connection.transaction()?;
    write_transaction(&transaction, system)?;
    transaction.commit()?;
    Ok(())
}

fn reserve_destination(path: &Path) -> Result<(), DmsError> {
    match OpenOptions::new().write(true).create_new(true).open(path) {
        Ok(file) => {
            drop(file);
            Ok(())
        }
        Err(error) if error.kind() == ErrorKind::AlreadyExists => Err(DmsError::DestinationExists),
        Err(error) => Err(DmsError::Io(error)),
    }
}

fn write_transaction(transaction: &Transaction<'_>, system: &DmsSystem) -> Result<(), DmsError> {
    transaction.execute_batch(PARTICLE_SCHEMA)?;
    write_version(transaction, system)?;
    write_particles(transaction, system)?;
    write_bonds(transaction, &system.topology.bonds)?;
    write_cell(transaction, system)?;
    Ok(())
}

fn write_version(transaction: &Transaction<'_>, system: &DmsSystem) -> Result<(), DmsError> {
    let Some(version) = system.version else {
        return Ok(());
    };
    transaction.execute_batch(
        "CREATE TABLE dms_version (major INTEGER NOT NULL, minor INTEGER NOT NULL)",
    )?;
    transaction.execute(
        "INSERT INTO dms_version (major, minor) VALUES (?1, ?2)",
        params![version.major, version.minor],
    )?;
    Ok(())
}

fn write_particles(transaction: &Transaction<'_>, system: &DmsSystem) -> Result<(), DmsError> {
    let mut statement = transaction.prepare(PARTICLE_INSERT)?;
    for (index, particle) in system.topology.particles.iter().enumerate() {
        let id = i64::try_from(index).map_err(|_| DmsError::InvalidParticleIds)?;
        let position = system.frame.positions[index];
        let velocity = system.frame.velocities[index];
        let (vx, vy, vz) = match velocity {
            Some(vector) => (Some(vector[0]), Some(vector[1]), Some(vector[2])),
            None => (None, None, None),
        };
        statement.execute(params![
            id,
            particle.atomic_number,
            particle.component,
            particle.nonbonded_type,
            position[0],
            position[1],
            position[2],
            vx,
            vy,
            vz,
            particle.mass,
            particle.charge,
            particle.residue_id,
            particle.residue_name.as_deref(),
            particle.chain.as_deref(),
            particle.segment.as_deref(),
            particle.name.as_deref(),
            particle.insertion.as_deref(),
            particle.formal_charge,
            particle.occupancy,
            particle.b_factor,
            particle.temperature_group,
            particle.energy_group,
            particle.ligand_group,
            particle.bias_group,
        ])?;
    }
    Ok(())
}

fn write_bonds(transaction: &Transaction<'_>, bonds: &[DmsBond]) -> Result<(), DmsError> {
    let mut ordered = bonds.to_vec();
    ordered.sort_unstable_by_key(|bond| (bond.p0, bond.p1));
    let mut statement =
        transaction.prepare("INSERT INTO bond (p0, p1, \"order\") VALUES (?1, ?2, ?3)")?;
    for bond in ordered {
        statement.execute(params![bond.p0, bond.p1, bond.order])?;
    }
    Ok(())
}

fn write_cell(transaction: &Transaction<'_>, system: &DmsSystem) -> Result<(), DmsError> {
    let Some(cell) = system.frame.cell else {
        return Ok(());
    };
    transaction.execute_batch(
        "CREATE TABLE global_cell (
          id INTEGER PRIMARY KEY,
          x REAL NOT NULL,
          y REAL NOT NULL,
          z REAL NOT NULL
        )",
    )?;
    let mut statement =
        transaction.prepare("INSERT INTO global_cell (id, x, y, z) VALUES (?1, ?2, ?3, ?4)")?;
    for (id, vector) in cell.vectors.into_iter().enumerate() {
        let id = i64::try_from(id).map_err(|_| DmsError::InvalidCell)?;
        statement.execute(params![id, vector[0], vector[1], vector[2]])?;
    }
    Ok(())
}

fn validate(system: &DmsSystem) -> Result<(), DmsError> {
    let count = system.topology.particles.len();
    if system.frame.positions.len() != count || system.frame.velocities.len() != count {
        return Err(DmsError::ParticleCountMismatch);
    }
    for (index, (particle, position)) in system
        .topology
        .particles
        .iter()
        .zip(&system.frame.positions)
        .enumerate()
    {
        validate_particle(index, particle, *position, system.frame.velocities[index])?;
    }
    validate_bonds(&system.topology.bonds, count)?;
    if let Some(cell) = system.frame.cell
        && cell
            .vectors
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(DmsError::InvalidCell);
    }
    Ok(())
}

fn validate_particle(
    index: usize,
    particle: &DmsParticle,
    position: [f64; 3],
    velocity: Option<[f64; 3]>,
) -> Result<(), DmsError> {
    let invalid_position = position.iter().any(|value| !value.is_finite());
    let invalid_velocity =
        velocity.is_some_and(|values| values.iter().any(|value| !value.is_finite()));
    let optional_values = [
        particle.mass,
        particle.charge,
        particle.formal_charge,
        particle.occupancy,
        particle.b_factor,
    ];
    if invalid_position
        || invalid_velocity
        || optional_values
            .into_iter()
            .flatten()
            .any(|value| !value.is_finite())
    {
        let id = i64::try_from(index).map_err(|_| DmsError::InvalidParticleIds)?;
        return Err(DmsError::InvalidParticle { id });
    }
    Ok(())
}

fn validate_bonds(bonds: &[DmsBond], particles: usize) -> Result<(), DmsError> {
    let mut seen = BTreeSet::new();
    for bond in bonds {
        let p0 = i64::from(bond.p0);
        let p1 = i64::from(bond.p1);
        if bond.p0 >= bond.p1
            || usize::try_from(bond.p1).map_or(true, |index| index >= particles)
            || !bond.order.is_finite()
            || !seen.insert((bond.p0, bond.p1))
        {
            return Err(DmsError::InvalidBond { p0, p1 });
        }
    }
    Ok(())
}
