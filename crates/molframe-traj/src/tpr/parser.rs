//! TPR header and top-level section decoding.

use super::constants::{DIMENSIONS, FORMAT_SIZE_FIELD, SUPPORTED_VERSIONS};
use super::error::TprError;
use super::model::{TprHeader, TprTopology};
use super::topology::parse_topology;
use super::xdr::Decoder;

const VERSION_PREFIX: &str = "VERSION";
const RELEASE_TAG: &str = "release";
const MODERN_GENERATION: i32 = 27;

/// Reads the structural topology from a portable GROMACS run-input file.
///
/// Simulation parameters, velocities and forces are deliberately not exposed:
/// they belong to the simulation engine rather than the structural topology.
///
/// # Errors
///
/// Returns an error for malformed input, unsupported format versions, missing
/// topology, or any versioned section that cannot be decoded without guessing.
pub fn parse_tpr(bytes: &[u8]) -> Result<TprTopology, TprError> {
    let mut decoder = Decoder::new(bytes);
    let internal = read_header(&mut decoder)?;
    if internal.modern {
        reject_beta_serializer(&decoder, internal.body_size)?;
        decoder.set_modern();
    }
    skip_state_prefix(&mut decoder, &internal)?;
    if !internal.has_topology {
        return Err(TprError::MissingTopology);
    }
    let mut topology = parse_topology(&mut decoder, &internal.public)?;
    topology.header = internal.public;
    Ok(topology)
}

#[derive(Debug)]
struct InternalHeader {
    public: TprHeader,
    temperature_groups: usize,
    has_topology: bool,
    has_box: bool,
    modern: bool,
    body_size: Option<i64>,
}

fn read_header(decoder: &mut Decoder<'_>) -> Result<InternalHeader, TprError> {
    let producer = decoder.string("producer version")?;
    if !producer.starts_with(VERSION_PREFIX) {
        return Err(TprError::InvalidMagic);
    }
    let precision = decoder.i32()?;
    decoder.set_precision(precision)?;
    let version = decoder.i32()?;
    if !SUPPORTED_VERSIONS.contains(&version) {
        return Err(TprError::UnsupportedVersion(version));
    }
    let early_tag = if (77..=79).contains(&version) {
        decoder.i32()?;
        Some(decoder.string("file tag")?)
    } else {
        None
    };
    let generation = if version >= 26 { decoder.i32()? } else { 0 };
    let tag = match early_tag {
        Some(tag) => tag,
        None if version >= 81 => decoder.string("file tag")?,
        None => RELEASE_TAG.to_owned(),
    };
    let atom_count = decoder.count("atom count")?;
    let temperature_groups = if version >= 28 {
        decoder.count("temperature group count")?
    } else {
        0
    };
    if version < 62 {
        decoder.i32()?;
        decoder.real()?;
    }
    if version >= 79 {
        decoder.i32()?;
    }
    decoder.real()?;
    let _has_input_record = decoder.i32()? != 0;
    let has_topology = decoder.i32()? != 0;
    let has_coordinates = decoder.i32()? != 0;
    let _has_velocities = decoder.i32()? != 0;
    let _has_forces = decoder.i32()? != 0;
    let has_box = decoder.i32()? != 0;
    let modern = version >= FORMAT_SIZE_FIELD && generation >= MODERN_GENERATION;
    let body_size = if modern { Some(decoder.i64()?) } else { None };
    let precision = u8::try_from(precision).map_err(|_| TprError::SizeOverflow)?;
    Ok(InternalHeader {
        public: TprHeader {
            producer,
            format_version: version,
            generation,
            precision,
            atom_count,
            tag,
            has_coordinates,
        },
        temperature_groups,
        has_topology,
        has_box,
        modern,
        body_size,
    })
}

fn reject_beta_serializer(decoder: &Decoder<'_>, body_size: Option<i64>) -> Result<(), TprError> {
    let Some(body_size) = body_size else {
        return Ok(());
    };
    if body_size < 0 {
        return Err(TprError::InvalidCount {
            field: "body size",
            value: body_size,
            offset: decoder.offset(),
        });
    }
    let body_size = usize::try_from(body_size).map_err(|_| TprError::SizeOverflow)?;
    if body_size.checked_mul(4) == Some(decoder.remaining()) {
        return Err(TprError::UnsupportedBetaSerializer);
    }
    Ok(())
}

fn skip_state_prefix(decoder: &mut Decoder<'_>, header: &InternalHeader) -> Result<(), TprError> {
    if header.has_box {
        decoder.skip_reals(DIMENSIONS * DIMENSIONS * 3)?;
    }
    if header.temperature_groups > 0 {
        if header.public.format_version < 69 {
            decoder.skip_reals(header.temperature_groups)?;
        }
        decoder.skip_reals(header.temperature_groups)?;
    }
    Ok(())
}
