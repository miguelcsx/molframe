use super::{HoomdXmlError, parse_hoomd_xml};

const XML: &str = "<?xml version=\"1.0\"?><hoomd_xml version=\"1.6\"><configuration time_step=\"42\" dimensions=\"3\" natoms=\"3\"><box lx=\"10\" ly=\"11\" lz=\"12\" xy=\"0.2\"/><position>0 0 0 1 0 0 0 1 0</position><image>0 0 0 1 0 0 0 -1 0</image><velocity>1 2 3 4 5 6 7 8 9</velocity><acceleration>0 0 1 0 1 0 1 0 0</acceleration><orientation>1 0 0 0 0 1 0 0 0 0 1 0</orientation><type>O H H</type><mass>16 1 1</mass><charge>-0.8 0.4 0.4</charge><bond>OH 0 1 OH 0 2</bond><angle>HOH 1 0 2</angle><custom>alpha beta</custom></configuration></hoomd_xml>";

#[test]
fn particle_arrays_cell_topology_and_extensions_are_retained() {
    let config = parse_hoomd_xml(XML).expect("valid HOOMD XML should parse");
    assert_eq!(config.step, Some(42));
    assert_eq!(config.positions.len(), 3);
    assert_eq!(config.types[0].as_ref(), "O");
    assert_eq!(config.images[1], [1, 0, 0]);
    assert!(
        config.velocities[2]
            .iter()
            .zip([7.0, 8.0, 9.0])
            .all(|(actual, expected)| (actual - expected).abs() < f64::EPSILON)
    );
    assert!(
        config.accelerations[0]
            .iter()
            .zip([0.0, 0.0, 1.0])
            .all(|(actual, expected)| (actual - expected).abs() < f64::EPSILON)
    );
    assert!(
        config.orientations[1]
            .iter()
            .zip([0.0, 1.0, 0.0, 0.0])
            .all(|(actual, expected)| (actual - expected).abs() < f64::EPSILON)
    );
    assert_eq!(config.bonds[1].particles, [0, 2]);
    assert_eq!(config.angles[0].particles, [1, 0, 2]);
    assert_eq!(config.extensions["custom"].as_ref(), "alpha beta");
    assert!((config.cell.expect("cell should exist").tilt[0] - 0.2).abs() < f64::EPSILON);
}

#[test]
fn mismatched_particle_arrays_are_rejected() {
    let broken = XML.replace("<type>O H H</type>", "<type>O H</type>");
    assert!(matches!(
        parse_hoomd_xml(&broken),
        Err(HoomdXmlError::LengthMismatch)
    ));
}

#[test]
fn out_of_range_interaction_indices_are_rejected() {
    let broken = XML.replace("OH 0 2", "OH 0 3");
    assert!(matches!(
        parse_hoomd_xml(&broken),
        Err(HoomdXmlError::InvalidParticle)
    ));
}
