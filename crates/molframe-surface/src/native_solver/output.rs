//! Streaming output sizing before a shared-memory reservation and decoding.
use super::{NativeAtoms, NativeSurfaceError, NativeSurfaceOptions, NativeSurfaceProvider};
use crate::{
    IndexedSurfaceMesh, MolecularSurfaceMesh, SurfaceOutputError, read_msms_surface,
    read_surf_surface,
};
use molframe_core::execution::{ExecutionContext, MemoryReservation};
use std::{
    fs::File,
    io::{BufReader, Read},
    path::{Path, PathBuf},
};

pub(super) fn read(
    directory: &Path,
    provider: NativeSurfaceProvider,
    atoms: NativeAtoms<'_>,
    options: NativeSurfaceOptions,
    context: &ExecutionContext,
) -> Result<(MolecularSurfaceMesh, MemoryReservation), NativeSurfaceError> {
    let files = match provider {
        NativeSurfaceProvider::Surf => vec![(directory.join("atoms.tri"), None)],
        NativeSurfaceProvider::Msms => component_files(directory, options.maximum_components)?,
    };
    let (vertices, faces, memory) = reserve_output(&files, options, context)?;
    let rows = atoms.atoms.rows;
    let combined_vertices = if files.len() > 1 { vertices } else { 0 };
    let combined_faces = if files.len() > 1 { faces } else { 0 };
    let mut positions = Vec::with_capacity(combined_vertices);
    let mut normals = Vec::with_capacity(combined_vertices);
    let mut triangles = Vec::with_capacity(combined_faces);
    let mut atom_rows = Vec::with_capacity(combined_vertices);
    let mut discarded_degenerate_faces = 0u64;
    for (vertex_file, face_file) in &files {
        let mut component = match face_file {
            None => read_surf_surface(
                BufReader::new(File::open(vertex_file)?),
                rows,
                options.output,
                context.cancellation(),
            )?,
            Some(face_file) => read_msms_surface(
                BufReader::new(File::open(vertex_file)?),
                BufReader::new(File::open(face_file)?),
                rows,
                options.output,
                context.cancellation(),
            )?,
        };
        if provider == NativeSurfaceProvider::Msms {
            super::provenance::validate_isolated_sphere(&mut component, atoms, context)?;
        }
        if files.len() == 1 {
            return Ok((component, memory));
        }
        discarded_degenerate_faces += component.discarded_degenerate_faces;
        let offset =
            u32::try_from(positions.len()).map_err(|_| SurfaceOutputError::Limit("vertices"))?;
        for face in &mut component.geometry.faces {
            for index in &mut face.0 {
                *index = index
                    .checked_add(offset)
                    .ok_or(SurfaceOutputError::Limit("vertices"))?;
            }
        }
        positions.append(&mut component.geometry.vertices);
        normals.append(&mut component.geometry.vertex_normals);
        triangles.append(&mut component.geometry.faces);
        atom_rows.extend_from_slice(&component.atom_rows);
    }
    context.cancellation().check()?;
    let geometry = IndexedSurfaceMesh::from_analytic_normals(
        positions,
        triangles,
        normals,
        context.cancellation(),
    )?;
    let mesh = MolecularSurfaceMesh {
        geometry,
        atom_rows: atom_rows.into_boxed_slice(),
        discarded_degenerate_faces,
    };
    Ok((mesh, memory))
}

fn reserve_output(
    files: &[(PathBuf, Option<PathBuf>)],
    options: NativeSurfaceOptions,
    context: &ExecutionContext,
) -> Result<(usize, usize, MemoryReservation), NativeSurfaceError> {
    let mut vertices = 0usize;
    let mut faces = 0usize;
    let mut remaining = options.output.bytes;
    for (vertex_file, face_file) in files {
        let bytes = std::fs::metadata(vertex_file)?
            .len()
            .checked_add(match face_file {
                Some(path) => std::fs::metadata(path)?.len(),
                None => 0,
            })
            .ok_or(SurfaceOutputError::Limit("bytes"))?;
        remaining = remaining
            .checked_sub(bytes)
            .ok_or(SurfaceOutputError::Limit("bytes"))?;
        let (new_vertices, new_faces) = if let Some(face_file) = face_file {
            (
                count(vertex_file, 9, bytes, context)?,
                count(face_file, 5, bytes, context)?,
            )
        } else {
            let triangles = count(vertex_file, 19, bytes, context)?;
            (
                triangles
                    .checked_mul(3)
                    .ok_or(SurfaceOutputError::Limit("vertices"))?,
                triangles,
            )
        };
        vertices = vertices
            .checked_add(new_vertices)
            .ok_or(SurfaceOutputError::Limit("vertices"))?;
        faces = faces
            .checked_add(new_faces)
            .ok_or(SurfaceOutputError::Limit("faces"))?;
    }
    if faces > options.output.faces
        || (files.iter().any(|(_, face)| face.is_some()) && vertices > options.output.vertices)
    {
        return Err(SurfaceOutputError::Limit("geometry").into());
    }
    // Unindexed Surf records can share vertices. The decoder enforces the
    // limit on distinct analytical attributes, not on repeated triangle corners.
    vertices = vertices.min(options.output.vertices);
    // Covers growing component columns, combined output and flat topology scratch.
    let bytes = vertices
        .checked_mul(128)
        .and_then(|size| {
            faces
                .checked_mul(128)
                .and_then(|face_bytes| size.checked_add(face_bytes))
        })
        .ok_or(SurfaceOutputError::Limit("memory"))?;
    let memory = context.try_reserve(bytes)?;
    Ok((vertices, faces, memory))
}

fn count(
    path: &Path,
    fields: usize,
    maximum: u64,
    context: &ExecutionContext,
) -> Result<usize, NativeSurfaceError> {
    let mut reader = File::open(path)?;
    if reader.metadata()?.len() > maximum {
        return Err(SurfaceOutputError::Limit("bytes").into());
    }
    let mut buffer = [0u8; 8192];
    let mut tokens = 0usize;
    let mut in_word = false;
    loop {
        context.cancellation().check()?;
        let length = reader.read(&mut buffer)?;
        if length == 0 {
            break;
        }
        for &byte in &buffer[..length] {
            if byte.is_ascii_whitespace() {
                in_word = false;
            } else if !in_word {
                tokens = tokens
                    .checked_add(1)
                    .ok_or(SurfaceOutputError::Limit("records"))?;
                in_word = true;
            }
        }
    }
    if !tokens.is_multiple_of(fields) {
        return Err(SurfaceOutputError::Invalid("incomplete record").into());
    }
    Ok(tokens / fields)
}

fn component_files(
    directory: &Path,
    maximum: usize,
) -> Result<Vec<(PathBuf, Option<PathBuf>)>, NativeSurfaceError> {
    let mut indices = Vec::new();
    for entry in std::fs::read_dir(directory)? {
        let entry = entry?;
        let name = entry.file_name();
        let Some(name) = name.to_str() else {
            continue;
        };
        if let Some(index) = name
            .strip_prefix("surface_")
            .and_then(|name| name.strip_suffix(".vert"))
        {
            let index = index
                .parse::<u32>()
                .map_err(|_| SurfaceOutputError::Invalid("invalid component filename"))?;
            indices.push(index);
        }
    }
    if indices.len() >= maximum {
        return Err(SurfaceOutputError::Limit("components").into());
    }
    indices.sort_unstable();
    let mut files = vec![(
        directory.join("surface.vert"),
        Some(directory.join("surface.face")),
    )];
    if !files[0].0.try_exists()? {
        return Err(SurfaceOutputError::Invalid(
            "MSMS completed without a triangulated main surface component",
        )
        .into());
    }
    for (expected, index) in (1u32..).zip(indices) {
        if expected != index {
            return Err(SurfaceOutputError::Invalid("missing surface component").into());
        }
        files.push((
            directory.join(format!("surface_{index}.vert")),
            Some(directory.join(format!("surface_{index}.face"))),
        ));
    }
    for (_, face_file) in &files {
        if let Some(face_file) = face_file
            && !face_file.try_exists()?
        {
            return Err(SurfaceOutputError::Invalid(
                "MSMS produced surface vertices without their complete face stream",
            )
            .into());
        }
    }
    Ok(files)
}
