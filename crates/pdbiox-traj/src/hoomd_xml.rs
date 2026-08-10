//! HOOMD-blue XML topology and initial-coordinate reader.

use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::collections::BTreeMap;

/// Triclinic HOOMD box parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HoomdBox {
    /// Box lengths.
    pub lengths: [f64; 3],
    /// `xy`, `xz`, and `yz` tilt factors.
    pub tilt: [f64; 3],
}

/// One named particle interaction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HoomdInteraction<const N: usize> {
    /// Interaction type name.
    pub kind: Box<str>,
    /// Zero-based particle indices.
    pub particles: [u32; N],
}

/// One HOOMD configuration with topology arrays.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HoomdConfiguration {
    /// Simulation step.
    pub step: Option<u64>,
    /// Spatial dimensionality, normally two or three.
    pub dimensions: Option<u8>,
    /// Declared particle count.
    pub particle_count: usize,
    /// Periodic cell.
    pub cell: Option<HoomdBox>,
    /// Particle coordinates.
    pub positions: Vec<[f64; 3]>,
    /// Periodic image counters for unwrapping coordinates.
    pub images: Vec<[i64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle accelerations.
    pub accelerations: Vec<[f64; 3]>,
    /// Particle orientations as HOOMD quaternions.
    pub orientations: Vec<[f64; 4]>,
    /// Particle type names.
    pub types: Vec<Box<str>>,
    /// Optional masses.
    pub masses: Vec<f64>,
    /// Optional charges.
    pub charges: Vec<f64>,
    /// Optional diameters.
    pub diameters: Vec<f64>,
    /// Optional integer body identifiers.
    pub bodies: Vec<i64>,
    /// Bonds.
    pub bonds: Vec<HoomdInteraction<2>>,
    /// Angles.
    pub angles: Vec<HoomdInteraction<3>>,
    /// Dihedrals.
    pub dihedrals: Vec<HoomdInteraction<4>>,
    /// Impropers.
    pub impropers: Vec<HoomdInteraction<4>>,
    /// Unknown per-particle or extension elements retained as text.
    pub extensions: BTreeMap<Box<str>, Box<str>>,
}

/// Malformed or inconsistent HOOMD XML.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum HoomdXmlError {
    /// XML tokenization failed.
    #[error("invalid HOOMD XML: {0}")]
    Xml(#[from] quick_xml::Error),
    /// Required configuration attributes or elements are absent.
    #[error("HOOMD XML has no valid configuration")]
    MissingConfiguration,
    /// A numeric value is malformed.
    #[error("invalid HOOMD XML numeric value")]
    InvalidNumber,
    /// An array length disagrees with the particle count.
    #[error("HOOMD XML particle array length mismatch")]
    LengthMismatch,
    /// An interaction names a particle outside the configuration.
    #[error("HOOMD XML interaction references an absent particle")]
    InvalidParticle,
}

/// Parses the first `<configuration>` in a HOOMD-blue XML document.
///
/// # Errors
///
/// Returns an error for malformed XML or numbers, missing configuration,
/// inconsistent particle arrays, and invalid interaction indices.
pub fn parse_hoomd_xml(text: &str) -> Result<HoomdConfiguration, HoomdXmlError> {
    let mut reader = Reader::from_str(text);
    reader.config_mut().trim_text(true);
    let mut output = None;
    let mut active: Option<Box<str>> = None;
    loop {
        match reader.read_event()? {
            Event::Start(element) if element.name().as_ref() == b"configuration" => {
                output = Some(configuration(&element)?);
            }
            Event::Empty(element) if element.name().as_ref() == b"box" => {
                let Some(config) = &mut output else { continue };
                config.cell = Some(parse_box(&element)?);
            }
            Event::Start(element) if output.is_some() => {
                active = Some(
                    String::from_utf8_lossy(element.name().as_ref())
                        .into_owned()
                        .into(),
                );
            }
            Event::Text(content) => {
                let (Some(config), Some(name)) = (&mut output, &active) else {
                    continue;
                };
                let decoded = content.decode().map_err(|_| HoomdXmlError::InvalidNumber)?;
                populate(config, name, &decoded)?;
            }
            Event::End(element) if element.name().as_ref() == b"configuration" => break,
            Event::End(_) => active = None,
            Event::Eof => break,
            _ => {}
        }
    }
    let Some(config) = output else {
        return Err(HoomdXmlError::MissingConfiguration);
    };
    validate(config)
}

fn configuration(element: &BytesStart<'_>) -> Result<HoomdConfiguration, HoomdXmlError> {
    let mut output = HoomdConfiguration::default();
    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(|_| HoomdXmlError::MissingConfiguration)?;
        let value = String::from_utf8_lossy(&attribute.value);
        match attribute.key.as_ref() {
            b"time_step" => {
                output.step = Some(value.parse().map_err(|_| HoomdXmlError::InvalidNumber)?);
            }
            b"dimensions" => {
                output.dimensions = Some(value.parse().map_err(|_| HoomdXmlError::InvalidNumber)?);
            }
            b"natoms" => {
                output.particle_count = value.parse().map_err(|_| HoomdXmlError::InvalidNumber)?;
            }
            _ => {}
        }
    }
    Ok(output)
}

fn parse_box(element: &BytesStart<'_>) -> Result<HoomdBox, HoomdXmlError> {
    let mut values = BTreeMap::new();
    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(|_| HoomdXmlError::InvalidNumber)?;
        let key = String::from_utf8_lossy(attribute.key.as_ref()).into_owned();
        let value = String::from_utf8_lossy(&attribute.value)
            .parse()
            .map_err(|_| HoomdXmlError::InvalidNumber)?;
        values.insert(key, value);
    }
    let required = |name: &str| {
        values
            .get(name)
            .copied()
            .ok_or(HoomdXmlError::InvalidNumber)
    };
    Ok(HoomdBox {
        lengths: [required("lx")?, required("ly")?, required("lz")?],
        tilt: [
            value_or_zero(&values, "xy"),
            value_or_zero(&values, "xz"),
            value_or_zero(&values, "yz"),
        ],
    })
}

fn value_or_zero(values: &BTreeMap<String, f64>, key: &str) -> f64 {
    match values.get(key) {
        Some(value) => *value,
        None => 0.0,
    }
}

fn populate(config: &mut HoomdConfiguration, name: &str, text: &str) -> Result<(), HoomdXmlError> {
    match name {
        "position" => config.positions = triples(text)?,
        "image" => config.images = integer_triples(text)?,
        "velocity" => config.velocities = triples(text)?,
        "acceleration" => config.accelerations = triples(text)?,
        "orientation" => config.orientations = quadruples(text)?,
        "type" => config.types = text.split_whitespace().map(Into::into).collect(),
        "mass" => config.masses = numbers(text)?,
        "charge" => config.charges = numbers(text)?,
        "diameter" => config.diameters = numbers(text)?,
        "body" => config.bodies = integers(text)?,
        "bond" => config.bonds = interactions(text)?,
        "angle" => config.angles = interactions(text)?,
        "dihedral" => config.dihedrals = interactions(text)?,
        "improper" => config.impropers = interactions(text)?,
        "box" => {}
        _ => {
            config.extensions.insert(name.into(), text.into());
        }
    }
    Ok(())
}

fn numbers(text: &str) -> Result<Vec<f64>, HoomdXmlError> {
    text.split_whitespace()
        .map(|value| value.parse().map_err(|_| HoomdXmlError::InvalidNumber))
        .collect()
}
fn integers(text: &str) -> Result<Vec<i64>, HoomdXmlError> {
    text.split_whitespace()
        .map(|value| value.parse().map_err(|_| HoomdXmlError::InvalidNumber))
        .collect()
}
fn triples(text: &str) -> Result<Vec<[f64; 3]>, HoomdXmlError> {
    let values = numbers(text)?;
    if values.len() % 3 != 0 {
        return Err(HoomdXmlError::LengthMismatch);
    }
    Ok(values
        .chunks_exact(3)
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

fn integer_triples(text: &str) -> Result<Vec<[i64; 3]>, HoomdXmlError> {
    let values = integers(text)?;
    if values.len() % 3 != 0 {
        return Err(HoomdXmlError::LengthMismatch);
    }
    Ok(values
        .chunks_exact(3)
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

fn quadruples(text: &str) -> Result<Vec<[f64; 4]>, HoomdXmlError> {
    let values = numbers(text)?;
    if values.len() % 4 != 0 {
        return Err(HoomdXmlError::LengthMismatch);
    }
    Ok(values
        .chunks_exact(4)
        .map(|row| [row[0], row[1], row[2], row[3]])
        .collect())
}

fn interactions<const N: usize>(text: &str) -> Result<Vec<HoomdInteraction<N>>, HoomdXmlError> {
    let fields: Vec<_> = text.split_whitespace().collect();
    if fields.len() % (N + 1) != 0 {
        return Err(HoomdXmlError::LengthMismatch);
    }
    fields
        .chunks_exact(N + 1)
        .map(|row| {
            let mut particles = [0; N];
            for (index, particle) in particles.iter_mut().enumerate() {
                *particle = row[index + 1]
                    .parse()
                    .map_err(|_| HoomdXmlError::InvalidNumber)?;
            }
            Ok(HoomdInteraction {
                kind: row[0].into(),
                particles,
            })
        })
        .collect()
}

fn validate(config: HoomdConfiguration) -> Result<HoomdConfiguration, HoomdXmlError> {
    if config.particle_count == 0
        || config.positions.len() != config.particle_count
        || config.types.len() != config.particle_count
    {
        return Err(HoomdXmlError::LengthMismatch);
    }
    for length in [
        config.images.len(),
        config.velocities.len(),
        config.accelerations.len(),
        config.orientations.len(),
        config.masses.len(),
        config.charges.len(),
        config.diameters.len(),
        config.bodies.len(),
    ] {
        if length != 0 && length != config.particle_count {
            return Err(HoomdXmlError::LengthMismatch);
        }
    }
    let invalid = config
        .bonds
        .iter()
        .any(|item| invalid_particles(&item.particles, config.particle_count))
        || config
            .angles
            .iter()
            .any(|item| invalid_particles(&item.particles, config.particle_count))
        || config
            .dihedrals
            .iter()
            .any(|item| invalid_particles(&item.particles, config.particle_count))
        || config
            .impropers
            .iter()
            .any(|item| invalid_particles(&item.particles, config.particle_count));
    if invalid {
        return Err(HoomdXmlError::InvalidParticle);
    }
    Ok(config)
}

fn invalid_particles<const N: usize>(particles: &[u32; N], count: usize) -> bool {
    particles
        .iter()
        .any(|particle| usize::try_from(*particle).map_or(true, |particle| particle >= count))
}

#[cfg(test)]
#[path = "hoomd_xml_tests.rs"]
mod tests;
