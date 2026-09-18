//! Bounded, exact IMD v2 wire encoding and decoding.

use std::io::{Read, Write};

use super::{ImdEnergies, ImdError, ImdForce, ImdLimits, ImdMessage, ImdPeerEndian};

pub(crate) const HEADER_BYTES: usize = 8;
pub(crate) const IMD_VERSION: i32 = 2;
const KILOCALORIE_TO_KILOJOULE: f32 = 4.184;
const ENERGY_SCALARS: usize = 9;

pub(crate) const DISCONNECT: i32 = 0;
pub(crate) const ENERGIES: i32 = 1;
pub(crate) const FCOORDS: i32 = 2;
pub(crate) const GO: i32 = 3;
pub(crate) const HANDSHAKE: i32 = 4;
pub(crate) const KILL: i32 = 5;
pub(crate) const MDCOMM: i32 = 6;
pub(crate) const PAUSE: i32 = 7;
pub(crate) const TRATE: i32 = 8;

pub(crate) fn negotiate(reader: &mut impl Read) -> Result<ImdPeerEndian, ImdError> {
    let mut header = [0u8; HEADER_BYTES];
    reader.read_exact(&mut header)?;
    if i32::from_be_bytes(first_word(header)) != HANDSHAKE {
        return Err(ImdError::InvalidHandshake);
    }
    let version = second_word(header);
    let same_endian = version == IMD_VERSION.to_ne_bytes();
    let opposite_endian = version == reversed(IMD_VERSION.to_ne_bytes());
    if !same_endian && !opposite_endian {
        return Err(ImdError::InvalidHandshake);
    }
    Ok(peer_endian(same_endian))
}

pub(crate) fn read_message(
    reader: &mut impl Read,
    endian: ImdPeerEndian,
    limits: ImdLimits,
) -> Result<ImdMessage, ImdError> {
    let (message_type, length) = read_header(reader)?;
    match message_type {
        DISCONNECT => control(length, message_type, ImdMessage::Disconnect),
        ENERGIES => read_energies(reader, endian, length),
        FCOORDS => read_coordinates(reader, endian, length, limits),
        GO => control(length, message_type, ImdMessage::Go),
        KILL => control(length, message_type, ImdMessage::Kill),
        MDCOMM => read_forces(reader, endian, length, limits),
        PAUSE => control(length, message_type, ImdMessage::Pause),
        TRATE => {
            let rate = u32::try_from(length).map_err(|_| invalid_length(message_type, length))?;
            if rate == 0 {
                Err(invalid_length(message_type, length))
            } else {
                Ok(ImdMessage::TransmissionRate(rate))
            }
        }
        _ => Err(ImdError::UnknownMessage(message_type)),
    }
}

pub(crate) fn write_control(
    writer: &mut impl Write,
    message_type: i32,
    length: i32,
) -> Result<(), ImdError> {
    writer.write_all(&header(message_type, length))?;
    writer.flush()?;
    Ok(())
}

pub(crate) fn write_forces(
    writer: &mut impl Write,
    endian: ImdPeerEndian,
    forces: &[ImdForce],
    limits: ImdLimits,
) -> Result<(), ImdError> {
    let count = i32::try_from(forces.len()).map_err(|_| ImdError::InvalidValue)?;
    if forces.len() > limits.max_atoms {
        return Err(invalid_length(MDCOMM, count));
    }
    validate_forces(forces)?;
    writer.write_all(&header(MDCOMM, count))?;
    for force in forces {
        let atom = i32::try_from(force.atom).map_err(|_| ImdError::InvalidValue)?;
        writer.write_all(&scalar_bytes(atom.cast_unsigned(), endian))?;
    }
    for force in forces {
        for component in force.force {
            let stored = component / KILOCALORIE_TO_KILOJOULE;
            writer.write_all(&scalar_bytes(stored.to_bits(), endian))?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn read_header(reader: &mut impl Read) -> Result<(i32, i32), ImdError> {
    let mut bytes = [0u8; HEADER_BYTES];
    reader.read_exact(&mut bytes)?;
    Ok((
        i32::from_be_bytes(first_word(bytes)),
        i32::from_be_bytes(second_word(bytes)),
    ))
}

fn read_energies(
    reader: &mut impl Read,
    endian: ImdPeerEndian,
    length: i32,
) -> Result<ImdMessage, ImdError> {
    if length != 1 {
        return Err(invalid_length(ENERGIES, length));
    }
    let step = read_i32(reader, endian)?;
    let temperature = read_f32(reader, endian)?;
    let mut energies = [0.0f32; ENERGY_SCALARS - 1];
    for energy in &mut energies {
        *energy = read_f32(reader, endian)? * KILOCALORIE_TO_KILOJOULE;
    }
    if !temperature.is_finite() || energies.iter().any(|value| !value.is_finite()) {
        return Err(ImdError::InvalidValue);
    }
    Ok(ImdMessage::Energies(ImdEnergies {
        step,
        temperature,
        total: energies[0],
        potential: energies[1],
        van_der_waals: energies[2],
        electrostatic: energies[3],
        bond: energies[4],
        angle: energies[5],
        dihedral: energies[6],
        improper: energies[7],
    }))
}

fn read_coordinates(
    reader: &mut impl Read,
    endian: ImdPeerEndian,
    length: i32,
    limits: ImdLimits,
) -> Result<ImdMessage, ImdError> {
    let count = checked_count(FCOORDS, length, limits)?;
    let mut positions = Vec::with_capacity(count);
    for _ in 0..count {
        let position = read_vector(reader, endian)?;
        if position.iter().any(|value| !value.is_finite()) {
            return Err(ImdError::InvalidValue);
        }
        positions.push(position);
    }
    Ok(ImdMessage::Coordinates(positions))
}

fn read_forces(
    reader: &mut impl Read,
    endian: ImdPeerEndian,
    length: i32,
    limits: ImdLimits,
) -> Result<ImdMessage, ImdError> {
    let count = checked_count(MDCOMM, length, limits)?;
    let mut atoms = Vec::with_capacity(count);
    for _ in 0..count {
        let atom = u32::try_from(read_i32(reader, endian)?).map_err(|_| ImdError::InvalidValue)?;
        atoms.push(atom);
    }
    let mut forces = Vec::with_capacity(count);
    for atom in atoms {
        let mut force = read_vector(reader, endian)?;
        for component in &mut force {
            *component *= KILOCALORIE_TO_KILOJOULE;
        }
        if force.iter().any(|value| !value.is_finite()) {
            return Err(ImdError::InvalidValue);
        }
        forces.push(ImdForce { atom, force });
    }
    Ok(ImdMessage::Forces(forces))
}

fn read_vector(reader: &mut impl Read, endian: ImdPeerEndian) -> Result<[f32; 3], ImdError> {
    Ok([
        read_f32(reader, endian)?,
        read_f32(reader, endian)?,
        read_f32(reader, endian)?,
    ])
}

fn read_i32(reader: &mut impl Read, endian: ImdPeerEndian) -> Result<i32, ImdError> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(match endian {
        ImdPeerEndian::Little => i32::from_le_bytes(bytes),
        ImdPeerEndian::Big => i32::from_be_bytes(bytes),
    })
}

fn read_f32(reader: &mut impl Read, endian: ImdPeerEndian) -> Result<f32, ImdError> {
    Ok(f32::from_bits(read_i32(reader, endian)?.cast_unsigned()))
}

fn checked_count(message_type: i32, length: i32, limits: ImdLimits) -> Result<usize, ImdError> {
    let count = usize::try_from(length).map_err(|_| invalid_length(message_type, length))?;
    if count > limits.max_atoms {
        Err(invalid_length(message_type, length))
    } else {
        Ok(count)
    }
}

fn control(length: i32, message_type: i32, message: ImdMessage) -> Result<ImdMessage, ImdError> {
    if length == 0 {
        Ok(message)
    } else {
        Err(invalid_length(message_type, length))
    }
}

fn validate_forces(forces: &[ImdForce]) -> Result<(), ImdError> {
    if forces.iter().all(|force| {
        i32::try_from(force.atom).is_ok() && force.force.iter().all(|value| value.is_finite())
    }) {
        Ok(())
    } else {
        Err(ImdError::InvalidValue)
    }
}

fn header(message_type: i32, length: i32) -> [u8; HEADER_BYTES] {
    let mut bytes = [0u8; HEADER_BYTES];
    bytes[..4].copy_from_slice(&message_type.to_be_bytes());
    bytes[4..].copy_from_slice(&length.to_be_bytes());
    bytes
}

fn scalar_bytes(bits: u32, endian: ImdPeerEndian) -> [u8; 4] {
    match endian {
        ImdPeerEndian::Little => bits.to_le_bytes(),
        ImdPeerEndian::Big => bits.to_be_bytes(),
    }
}

fn first_word(bytes: [u8; HEADER_BYTES]) -> [u8; 4] {
    [bytes[0], bytes[1], bytes[2], bytes[3]]
}

fn second_word(bytes: [u8; HEADER_BYTES]) -> [u8; 4] {
    [bytes[4], bytes[5], bytes[6], bytes[7]]
}

fn reversed(mut bytes: [u8; 4]) -> [u8; 4] {
    bytes.reverse();
    bytes
}

const fn peer_endian(same_endian: bool) -> ImdPeerEndian {
    match (cfg!(target_endian = "little"), same_endian) {
        (true, true) | (false, false) => ImdPeerEndian::Little,
        (false, true) | (true, false) => ImdPeerEndian::Big,
    }
}

const fn invalid_length(message_type: i32, length: i32) -> ImdError {
    ImdError::InvalidLength {
        message_type,
        length,
    }
}
