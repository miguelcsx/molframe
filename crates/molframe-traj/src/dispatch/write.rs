//! Path-based trajectory writing.

use std::path::Path;

use crate::numeric::f32_from_f64;

use super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryWriteOptions,
};

/// Writes normalized trajectory data by suffix or explicit format override.
///
/// # Errors
///
/// Returns an error when required steps/units are absent, metadata cannot be
/// represented, filesystem access fails, or the selected encoder refuses data.
pub fn write_trajectory(
    path: &Path,
    trajectory: &TrajectoryData,
    options: &TrajectoryWriteOptions,
) -> Result<(), TrajectoryIoError> {
    let format = match options.format.or_else(|| TrajectoryFormat::infer(path)) {
        Some(format) => format,
        None => trajectory.format,
    };
    let frames = frames_with_steps(trajectory, format)?;
    match format {
        TrajectoryFormat::Xtc => {
            let bytes = match options.xtc {
                Some(writer) => crate::write_xtc(&frames, writer)?,
                None => match xtc_precisions(&trajectory.metadata.format, frames.len())? {
                    Some(precisions) => crate::write_xtc_with_precisions(&frames, precisions)?,
                    None => crate::write_xtc(&frames, crate::XtcWriteOptions::default())?,
                },
            };
            std::fs::write(path, bytes)?;
        }
        TrajectoryFormat::Trr => {
            let bytes = match options.trr {
                Some(writer) => crate::write_trr(&frames, writer)?,
                None => match trr_precisions(&trajectory.metadata.format, frames.len())? {
                    Some(precisions) => crate::write_trr_with_precisions(&frames, precisions)?,
                    None => crate::write_trr(&frames, crate::TrrWriteOptions::default())?,
                },
            };
            std::fs::write(path, bytes)?;
        }
        TrajectoryFormat::Dcd => {
            let writer = match &options.dcd {
                Some(writer) => writer.clone(),
                None => dcd_options(&trajectory.metadata.format)?,
            };
            std::fs::write(path, crate::write_dcd(&frames, &writer)?)?;
        }
        TrajectoryFormat::AmberNetcdf => {
            let writer = amber_options_or_default(
                options
                    .amber_netcdf
                    .clone()
                    .or(amber_netcdf_options(&trajectory.metadata.format)?),
            );
            std::fs::write(path, crate::write_amber_netcdf(&frames, writer)?)?;
        }
        TrajectoryFormat::Tng => {
            let writer = match options.tng {
                Some(writer) => writer,
                None => tng_options(&trajectory.metadata.format),
            };
            crate::write_tng(path, &frames, writer)?;
        }
        TrajectoryFormat::H5md => {
            let bytes = match (&options.h5md, &trajectory.metadata.format) {
                (Some(h5md), _) => crate::write_h5md(&frames, h5md)?,
                (None, FormatMetadata::H5md(metadata)) => {
                    crate::write_h5md_with_metadata(&frames, metadata)?
                }
                (None, _) => crate::write_h5md(&frames, &crate::H5mdOptions::default())?,
            };
            std::fs::write(path, bytes)?;
        }
        TrajectoryFormat::Trz => {
            let title = match &options.trz {
                Some(options) => options.title.clone(),
                None => match &trajectory.metadata.format {
                    FormatMetadata::Trz { title, .. } => title.clone(),
                    _ => return Err(TrajectoryIoError::InvalidMetadata),
                },
            };
            let has_forces = frames.iter().all(|frame| frame.forces.is_some());
            let source = crate::TrzTrajectory {
                title,
                has_forces,
                frames,
            };
            std::fs::write(path, crate::write_trz(&source)?)?;
        }
        TrajectoryFormat::Namd | TrajectoryFormat::AmberRestart => {
            write_snapshot(path, trajectory, options, &frames, format)?;
        }
        TrajectoryFormat::Gro => write_gro_path(path, trajectory, &frames)?,
        TrajectoryFormat::Xyz => write_xyz_path(path, trajectory, &frames)?,
        TrajectoryFormat::Aims => write_aims_path(path, trajectory, &frames)?,
        TrajectoryFormat::Txyz => write_txyz_path(path, trajectory, &frames)?,
        TrajectoryFormat::DlPolyConfig | TrajectoryFormat::DlPolyHistory => {
            write_dlpoly_path(path, trajectory, &frames, format)?;
        }
        TrajectoryFormat::CharmmCard => write_charmm_path(path, trajectory, &frames)?,
        TrajectoryFormat::AmberAscii
        | TrajectoryFormat::Gsd
        | TrajectoryFormat::Gamess
        | TrajectoryFormat::LammpsDump
        | TrajectoryFormat::Gromos11
        | TrajectoryFormat::Dms => {
            return Err(TrajectoryIoError::ReadOnlyFormat);
        }
    }
    Ok(())
}

fn amber_options_or_default(
    options: Option<crate::AmberNetcdfWriteOptions>,
) -> crate::AmberNetcdfWriteOptions {
    let Some(options) = options else {
        return crate::AmberNetcdfWriteOptions::default();
    };
    options
}

fn amber_netcdf_options(
    metadata: &FormatMetadata,
) -> Result<Option<crate::AmberNetcdfWriteOptions>, TrajectoryIoError> {
    let FormatMetadata::AmberNetcdf(source) = metadata else {
        return Ok(None);
    };
    Ok(Some(crate::AmberNetcdfWriteOptions {
        precision: source.precision,
        program: source
            .program
            .as_deref()
            .ok_or(TrajectoryIoError::InvalidMetadata)?
            .into(),
        program_version: source
            .program_version
            .as_deref()
            .ok_or(TrajectoryIoError::InvalidMetadata)?
            .into(),
    }))
}

fn write_charmm_path(
    path: &Path,
    trajectory: &TrajectoryData,
    frames: &[crate::Timestep],
) -> Result<(), TrajectoryIoError> {
    let ([frame], FormatMetadata::CharmmCard(source)) = (frames, &trajectory.metadata.format)
    else {
        return Err(TrajectoryIoError::InvalidMetadata);
    };
    if source.atoms.len() != frame.positions.len()
        || frame.velocities.is_some()
        || frame.forces.is_some()
        || frame.cell.is_some()
        || frame.time.is_some()
        || frame.dt.is_some()
        || !frame.data.is_empty()
    {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    let mut card = source.clone();
    for (atom, position) in card.atoms.iter_mut().zip(&frame.positions) {
        atom.position = *position;
    }
    std::fs::write(path, crate::write_charmm_card(&card)?)?;
    Ok(())
}

fn write_dlpoly_path(
    path: &Path,
    trajectory: &TrajectoryData,
    frames: &[crate::Timestep],
    format: TrajectoryFormat,
) -> Result<(), TrajectoryIoError> {
    let rendered = match (format, &trajectory.metadata.format) {
        (TrajectoryFormat::DlPolyConfig, FormatMetadata::DlPolyConfig(source)) => {
            let [frame] = frames else {
                return Err(TrajectoryIoError::InvalidMetadata);
            };
            let mut source = source.clone();
            sync_dlpoly_frame(&mut source.frame, frame)?;
            crate::write_dlpoly_config(&source)?
        }
        (TrajectoryFormat::DlPolyHistory, FormatMetadata::DlPolyHistory(source)) => {
            if source.frames.len() != frames.len() {
                return Err(TrajectoryIoError::InvalidMetadata);
            }
            let mut source = source.clone();
            for (target, frame) in source.frames.iter_mut().zip(frames) {
                sync_dlpoly_frame(target, frame)?;
            }
            crate::write_dlpoly_history(&source)?
        }
        _ => return Err(TrajectoryIoError::InvalidMetadata),
    };
    std::fs::write(path, rendered)?;
    Ok(())
}

fn sync_dlpoly_frame(
    target: &mut crate::DlPolyFrame,
    frame: &crate::Timestep,
) -> Result<(), TrajectoryIoError> {
    if frame.cell
        != target
            .lattice_vectors
            .and_then(crate::cell::cell_from_vectors)
        || frame.dt.is_some()
        || frame.data.len() != 1
        || frame.data.get("dlpoly.step") != Some(&crate::FrameValue::Integer(target.step))
    {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    target.positions.clone_from(&frame.positions);
    target.velocities.clone_from(&frame.velocities);
    target.forces.clone_from(&frame.forces);
    target.time = frame.time.ok_or(TrajectoryIoError::InvalidMetadata)?;
    Ok(())
}

fn write_txyz_path(
    path: &Path,
    trajectory: &TrajectoryData,
    frames: &[crate::Timestep],
) -> Result<(), TrajectoryIoError> {
    let FormatMetadata::Txyz(source) = &trajectory.metadata.format else {
        return Err(TrajectoryIoError::InvalidMetadata);
    };
    if source.len() != frames.len() {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    let mut records = source.clone();
    for (record, frame) in records.iter_mut().zip(frames) {
        if record.atoms.len() != frame.positions.len()
            || frame.velocities.is_some()
            || frame.forces.is_some()
            || frame.cell.is_some()
            || frame.time.is_some()
            || frame.dt.is_some()
            || !frame.data.is_empty()
        {
            return Err(TrajectoryIoError::InvalidMetadata);
        }
        for (atom, position) in record.atoms.iter_mut().zip(&frame.positions) {
            atom.position = *position;
        }
    }
    std::fs::write(path, crate::write_txyz(&records)?)?;
    Ok(())
}

fn write_aims_path(
    path: &Path,
    trajectory: &TrajectoryData,
    frames: &[crate::Timestep],
) -> Result<(), TrajectoryIoError> {
    let ([frame], FormatMetadata::Aims(source)) = (frames, &trajectory.metadata.format) else {
        return Err(TrajectoryIoError::InvalidMetadata);
    };
    if source.atoms.len() != frame.positions.len()
        || frame.cell != source.cell()
        || frame.velocities.is_some()
        || frame.forces.is_some()
        || frame.time.is_some()
        || frame.dt.is_some()
        || !frame.data.is_empty()
    {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    let mut geometry = source.clone();
    for (atom, position) in geometry.atoms.iter_mut().zip(&frame.positions) {
        atom.position = *position;
    }
    std::fs::write(path, crate::write_aims_geometry(&geometry))?;
    Ok(())
}

fn write_xyz_path(
    path: &Path,
    trajectory: &TrajectoryData,
    frames: &[crate::Timestep],
) -> Result<(), TrajectoryIoError> {
    let FormatMetadata::Xyz(source) = &trajectory.metadata.format else {
        return Err(TrajectoryIoError::InvalidMetadata);
    };
    if source.len() != frames.len() {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    let mut records = source.clone();
    for (record, frame) in records.iter_mut().zip(frames) {
        if record.atoms.len() != frame.positions.len()
            || frame.velocities.is_some()
            || frame.forces.is_some()
            || frame.cell.is_some()
            || frame.time.is_some()
            || frame.dt.is_some()
            || !frame.data.is_empty()
        {
            return Err(TrajectoryIoError::InvalidMetadata);
        }
        for (atom, position) in record.atoms.iter_mut().zip(&frame.positions) {
            atom.position = *position;
        }
    }
    std::fs::write(path, crate::write_xyz(&records))?;
    Ok(())
}

fn write_gro_path(
    path: &Path,
    trajectory: &TrajectoryData,
    frames: &[crate::Timestep],
) -> Result<(), TrajectoryIoError> {
    let FormatMetadata::Gro(source) = &trajectory.metadata.format else {
        return Err(TrajectoryIoError::InvalidMetadata);
    };
    if source.len() != frames.len() {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    let mut records = source.clone();
    for (record, frame) in records.iter_mut().zip(frames) {
        if record.atoms.len() != frame.positions.len()
            || frame.forces.is_some()
            || frame.time.is_some()
            || frame.dt.is_some()
            || !frame.data.is_empty()
        {
            return Err(TrajectoryIoError::InvalidMetadata);
        }
        let velocities = match &frame.velocities {
            Some(values) if values.len() == record.atoms.len() => Some(values),
            None => None,
            _ => return Err(TrajectoryIoError::InvalidMetadata),
        };
        for (index, atom) in record.atoms.iter_mut().enumerate() {
            atom.position = frame.positions[index];
            atom.velocity = velocities.map(|values| values[index]);
        }
    }
    std::fs::write(path, crate::write_gro(&records)?)?;
    Ok(())
}

fn write_snapshot(
    path: &Path,
    trajectory: &TrajectoryData,
    options: &TrajectoryWriteOptions,
    frames: &[crate::Timestep],
    format: TrajectoryFormat,
) -> Result<(), TrajectoryIoError> {
    let [frame] = frames else {
        return Err(TrajectoryIoError::InvalidMetadata);
    };
    let bytes = match format {
        TrajectoryFormat::Namd => {
            let endian = options
                .namd
                .or(match trajectory.metadata.format {
                    FormatMetadata::Namd(endian) => Some(endian),
                    _ => None,
                })
                .ok_or(TrajectoryIoError::InvalidMetadata)?;
            crate::write_namd_binary(frame, endian)?
        }
        TrajectoryFormat::AmberRestart => {
            let FormatMetadata::AmberRestart { title, layout } = &trajectory.metadata.format else {
                return Err(TrajectoryIoError::InvalidMetadata);
            };
            crate::write_amber_restart(&crate::AmberRestart {
                title: title.clone(),
                layout: *layout,
                timestep: frame.clone(),
            })?
            .into_bytes()
        }
        _ => return Err(TrajectoryIoError::UnknownFormat),
    };
    std::fs::write(path, bytes)?;
    Ok(())
}

fn frames_with_steps(
    trajectory: &TrajectoryData,
    format: TrajectoryFormat,
) -> Result<Vec<crate::Timestep>, TrajectoryIoError> {
    if !matches!(
        format,
        TrajectoryFormat::Xtc
            | TrajectoryFormat::Trr
            | TrajectoryFormat::Tng
            | TrajectoryFormat::Gsd
            | TrajectoryFormat::H5md
    ) {
        return Ok(trajectory.frames.clone());
    }
    let steps = trajectory
        .metadata
        .steps
        .as_ref()
        .ok_or(TrajectoryIoError::InvalidMetadata)?;
    if steps.len() != trajectory.frames.len() {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    trajectory
        .frames
        .iter()
        .cloned()
        .zip(steps)
        .map(|(mut frame, step)| {
            frame.frame = usize::try_from(*step).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
            Ok(frame)
        })
        .collect()
}

fn xtc_precisions(
    metadata: &FormatMetadata,
    frames: usize,
) -> Result<Option<&[f32]>, TrajectoryIoError> {
    let FormatMetadata::Xtc { precision } = metadata else {
        return Ok(None);
    };
    if precision.len() != frames {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    Ok(Some(precision))
}

fn trr_precisions(
    metadata: &FormatMetadata,
    frames: usize,
) -> Result<Option<&[crate::TrrPrecision]>, TrajectoryIoError> {
    let FormatMetadata::Trr { precision } = metadata else {
        return Ok(None);
    };
    if precision.len() != frames {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    Ok(Some(precision))
}

fn dcd_options(metadata: &FormatMetadata) -> Result<crate::DcdWriteOptions, TrajectoryIoError> {
    let FormatMetadata::Dcd(header) = metadata else {
        return Ok(crate::DcdWriteOptions::default());
    };
    if !header.charmm
        || header.fixed_atom_count != 0
        || header.has_fourth_dimension
        || header.titles.len() != 1
    {
        return Err(TrajectoryIoError::InvalidMetadata);
    }
    let title = header
        .titles
        .first()
        .cloned()
        .ok_or(TrajectoryIoError::InvalidMetadata)?;
    Ok(crate::DcdWriteOptions {
        endian: header.endian,
        title,
        start_step: header.start_step,
        save_interval: header.save_interval,
        delta_akma: f32_from_f64(header.delta_akma)
            .ok_or(TrajectoryIoError::UnrepresentableValue)?,
    })
}

fn tng_options(metadata: &FormatMetadata) -> crate::TngWriteOptions {
    let FormatMetadata::Tng {
        distance_unit_exponent,
        compression,
        ..
    } = metadata
    else {
        return crate::TngWriteOptions::default();
    };
    crate::TngWriteOptions {
        distance_unit_exponent: *distance_unit_exponent,
        compression: *compression,
        ..crate::TngWriteOptions::default()
    }
}
