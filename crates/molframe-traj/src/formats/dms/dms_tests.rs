use super::*;
use rusqlite::Connection;

fn sample() -> DmsSystem {
    DmsSystem {
        version: Some(DmsVersion { major: 1, minor: 7 }),
        topology: DmsTopology {
            particles: vec![
                DmsParticle {
                    atomic_number: Some(6),
                    component: Some(2),
                    nonbonded_type: Some(4),
                    mass: Some(12.011),
                    charge: Some(-0.2),
                    residue_id: Some(42),
                    residue_name: Some("LIG".into()),
                    chain: Some("A".into()),
                    segment: Some("SOLUTE".into()),
                    name: Some("C1".into()),
                    insertion: Some("B".into()),
                    formal_charge: Some(0.0),
                    occupancy: Some(0.75),
                    b_factor: Some(12.5),
                    temperature_group: Some(1),
                    energy_group: Some(2),
                    ligand_group: Some(3),
                    bias_group: Some(4),
                },
                DmsParticle {
                    atomic_number: Some(8),
                    name: Some("O1".into()),
                    mass: Some(15.999),
                    charge: Some(-0.8),
                    ..DmsParticle::default()
                },
            ],
            bonds: vec![DmsBond {
                p0: 0,
                p1: 1,
                order: 1.5,
            }],
        },
        frame: DmsFrame {
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            velocities: vec![Some([0.1, 0.2, 0.3]), None],
            cell: Some(DmsCell {
                vectors: [[20.0, 0.0, 0.0], [2.0, 21.0, 0.0], [1.0, 3.0, 22.0]],
            }),
        },
    }
}

#[test]
fn writer_and_reader_round_trip_every_standard_structure_field() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("system.dms");
    let expected = sample();
    write_dms(&path, &expected).expect("write DMS");
    let actual = read_dms(&path).expect("read DMS");
    assert_eq!(actual, expected);
}

#[test]
fn minimal_schema_preserves_absent_optional_values_as_absent() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("minimal.dms");
    let connection = Connection::open(&path).expect("open SQLite");
    connection
        .execute_batch(
            "CREATE TABLE particle (
                id INTEGER PRIMARY KEY, x FLOAT, y FLOAT, z FLOAT
             );
             CREATE TABLE bond (p0 INTEGER, p1 INTEGER, 'order' FLOAT);
             INSERT INTO particle VALUES (0, 1.0, 2.0, 3.0);",
        )
        .expect("create minimal DMS");
    drop(connection);
    let system = read_dms(&path).expect("read minimal DMS");
    assert_eq!(system.topology.particles, vec![DmsParticle::default()]);
    assert_eq!(system.frame.positions, vec![[1.0, 2.0, 3.0]]);
    assert_eq!(system.frame.velocities, vec![None]);
    assert_eq!(system.frame.cell, None);
}

#[test]
fn non_contiguous_particle_identifiers_are_rejected() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("ids.dms");
    let connection = Connection::open(&path).expect("open SQLite");
    connection
        .execute_batch(
            "CREATE TABLE particle (
                id INTEGER PRIMARY KEY, x FLOAT, y FLOAT, z FLOAT
             );
             CREATE TABLE bond (p0 INTEGER, p1 INTEGER, 'order' FLOAT);
             INSERT INTO particle VALUES (1, 1.0, 2.0, 3.0);",
        )
        .expect("create malformed DMS");
    drop(connection);
    assert!(matches!(read_dms(&path), Err(DmsError::InvalidParticleIds)));
}

#[test]
fn partial_velocity_vectors_are_rejected() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("velocity.dms");
    let connection = Connection::open(&path).expect("open SQLite");
    connection
        .execute_batch(
            "CREATE TABLE particle (
                id INTEGER PRIMARY KEY, x FLOAT, y FLOAT, z FLOAT,
                vx FLOAT, vy FLOAT, vz FLOAT
             );
             CREATE TABLE bond (p0 INTEGER, p1 INTEGER, 'order' FLOAT);
             INSERT INTO particle VALUES (0, 1.0, 2.0, 3.0, 0.1, NULL, 0.3);",
        )
        .expect("create malformed DMS");
    drop(connection);
    assert!(matches!(
        read_dms(&path),
        Err(DmsError::InvalidParticle { id: 0 })
    ));
}

#[test]
fn invalid_or_duplicate_connectivity_is_rejected_before_writing() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("bonds.dms");
    let mut system = sample();
    system.topology.bonds.push(system.topology.bonds[0]);
    assert!(matches!(
        write_dms(&path, &system),
        Err(DmsError::InvalidBond { p0: 0, p1: 1 })
    ));
    assert!(!path.exists());
}

#[test]
fn writer_refuses_to_replace_an_existing_database() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("existing.dms");
    Connection::open(&path).expect("create destination");
    assert!(matches!(
        write_dms(&path, &sample()),
        Err(DmsError::DestinationExists)
    ));
}

#[test]
fn a_present_cell_requires_exact_zero_based_vector_ids() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("cell.dms");
    let connection = Connection::open(&path).expect("open SQLite");
    connection
        .execute_batch(
            "CREATE TABLE particle (
                id INTEGER PRIMARY KEY, x FLOAT, y FLOAT, z FLOAT
             );
             CREATE TABLE bond (p0 INTEGER, p1 INTEGER, 'order' FLOAT);
             CREATE TABLE global_cell (
                id INTEGER PRIMARY KEY, x FLOAT, y FLOAT, z FLOAT
             );
             INSERT INTO particle VALUES (0, 1.0, 2.0, 3.0);
             INSERT INTO global_cell VALUES (1, 1.0, 0.0, 0.0);
             INSERT INTO global_cell VALUES (2, 0.0, 1.0, 0.0);
             INSERT INTO global_cell VALUES (3, 0.0, 0.0, 1.0);",
        )
        .expect("create malformed DMS");
    drop(connection);
    assert!(matches!(read_dms(&path), Err(DmsError::InvalidCell)));
}

#[test]
fn projection_only_synthesizes_absent_optional_columns() {
    let columns = crate::dms::reader::column_set(&["id", "x"]);
    assert_eq!(
        crate::dms::schema::select_projection(&columns, &["id", "anum", "x"]),
        "id, NULL AS anum, x"
    );
}

#[test]
fn incompatible_declared_column_affinity_is_rejected() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("types.dms");
    let connection = Connection::open(&path).expect("open SQLite");
    connection
        .execute_batch(
            "CREATE TABLE particle (
                id TEXT PRIMARY KEY, x FLOAT, y FLOAT, z FLOAT
             );
             CREATE TABLE bond (p0 INTEGER, p1 INTEGER, 'order' FLOAT);
             INSERT INTO particle VALUES ('0', 1.0, 2.0, 3.0);",
        )
        .expect("create malformed DMS");
    drop(connection);
    assert!(matches!(
        read_dms(&path),
        Err(DmsError::InvalidColumnType {
            table: "particle",
            column: "id"
        })
    ));
}
