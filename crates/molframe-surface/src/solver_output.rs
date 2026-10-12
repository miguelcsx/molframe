//! Bounded decoding of native molecular-surface solver outputs.
//!
//! Parsing costs O(vertices + faces); shared mesh diagnostics cost O(faces log faces).
//! One reusable line buffer serves all records. Source atom identities and the
//! solver's analytic normals survive decoding rather than being guessed later.

use crate::{IndexedSurfaceMesh, SurfaceFace};
use molframe_core::execution::CancellationToken;
use num_traits::ToPrimitive;
use std::io::{self, BufRead, Read};

/// An indexed molecular surface retaining its provider's atom provenance.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularSurfaceMesh {
    /// Validated triangle geometry and the provider's unit vertex normals.
    pub geometry: IndexedSurfaceMesh,
    /// One original molecular atom row per vertex, in geometry order.
    pub atom_rows: Box<[u32]>,
    /// Provider triangles with exactly zero area, omitted without moving vertices.
    pub discarded_degenerate_faces: u64,
}

/// Hard limits on accepted surface output, independent of solver behavior.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SurfaceOutputLimits {
    /// Maximum decoded vertices.
    pub vertices: usize,
    /// Maximum decoded faces.
    pub faces: usize,
    /// Maximum textual bytes across all output files.
    pub bytes: u64,
}

impl Default for SurfaceOutputLimits {
    fn default() -> Self {
        Self {
            vertices: 4_000_000,
            faces: 8_000_000,
            bytes: 1_000_000_000,
        }
    }
}

/// A provider output could not be accepted as complete molecular geometry.
#[derive(Debug, thiserror::Error)]
pub enum SurfaceOutputError {
    /// The filesystem or input reader failed.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Cooperative cancellation was requested.
    #[error(transparent)]
    Cancelled(#[from] molframe_core::execution::Cancelled),
    /// The provider record shape, number, normal or atom identity is invalid.
    #[error("invalid molecular surface output: {0}")]
    Invalid(&'static str),
    /// A complete output cannot fit the caller's limits.
    #[error("molecular surface output exceeds its {0} limit")]
    Limit(&'static str),
}

/// Reads headerless MSMS vertex and face streams with one-based provider indices.
/// `input_atom_rows[i]` is the original atom row of input sphere `i`.
///
/// # Errors
/// Rejects malformed, nonfinite, out-of-range, truncated or over-budget output.
/// No incomplete prefix is returned.
pub fn read_msms_surface(
    vertices: impl BufRead,
    faces: impl BufRead,
    input_atom_rows: &[u32],
    limits: SurfaceOutputLimits,
    cancellation: &CancellationToken,
) -> Result<MolecularSurfaceMesh, SurfaceOutputError> {
    let mut vertex_reader = Records::new(vertices, limits.bytes, cancellation);
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut atom_rows = Vec::new();
    while let Some(record) = vertex_reader.next::<9>()? {
        if positions.len() >= limits.vertices {
            return Err(SurfaceOutputError::Limit("vertices"));
        }
        positions.push(point(&record[..3])?);
        normals.push(normal(&record[3..6])?);
        if record[6].fract() != 0.0 || record[6].to_i32().is_none() {
            return Err(SurfaceOutputError::Invalid(
                "invalid signed analytical face index",
            ));
        }
        let sphere = index(record[7])?
            .checked_sub(1)
            .ok_or(SurfaceOutputError::Invalid(
                "MSMS atom indices are one based",
            ))?;
        atom_rows.push(atom_row(input_atom_rows, sphere)?);
        let _ = index(record[8])?;
    }
    let mut face_reader = Records::new(
        faces,
        limits.bytes.saturating_sub(vertex_reader.bytes),
        cancellation,
    );
    let mut triangles = Vec::new();
    let mut decoded_faces = 0usize;
    let mut discarded = 0u64;
    while let Some(record) = face_reader.next::<5>()? {
        if decoded_faces >= limits.faces {
            return Err(SurfaceOutputError::Limit("faces"));
        }
        decoded_faces += 1;
        let mut face = [0u32; 3];
        for (out, &value) in face.iter_mut().zip(&record[..3]) {
            *out = index(value)?
                .checked_sub(1)
                .ok_or(SurfaceOutputError::Invalid(
                    "MSMS face indices are one based",
                ))?;
        }
        let _ = index(record[3])?;
        let _ = index(record[4])?;
        let [a, b, c] = face.map(|row| positions.get(row as usize));
        let (Some(a), Some(b), Some(c)) = (a, b, c) else {
            return Err(SurfaceOutputError::Invalid("face addresses no vertex"));
        };
        if crate::triangulation::triangle_cross(*a, *b, *c)
            .iter()
            .all(|&value| value == 0.0)
        {
            discarded += 1;
        } else {
            triangles.push(SurfaceFace(face));
        }
    }
    finish(
        positions,
        normals,
        atom_rows,
        triangles,
        discarded,
        cancellation,
    )
}

/// Reads Surf triangle records: one zero-based atom row followed by three
/// coordinate/normal pairs. Records may span physical lines.
///
/// # Errors
/// Rejects invalid or incomplete records and geometry beyond explicit limits.
pub fn read_surf_surface(
    input: impl BufRead,
    input_atom_rows: &[u32],
    limits: SurfaceOutputLimits,
    cancellation: &CancellationToken,
) -> Result<MolecularSurfaceMesh, SurfaceOutputError> {
    let mut reader = Records::new(input, limits.bytes, cancellation);
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut atom_rows = Vec::new();
    let mut faces = Vec::new();
    let mut vertices = std::collections::HashMap::new();
    let mut decoded_faces = 0usize;
    let mut discarded = 0u64;
    while let Some(record) = reader.next::<19>()? {
        if decoded_faces >= limits.faces {
            return Err(SurfaceOutputError::Limit("faces"));
        }
        decoded_faces += 1;
        let atom = atom_row(input_atom_rows, index(record[0])?)?;
        let mut triangle = [[0.0; 3]; 3];
        let mut directions = [[0.0; 3]; 3];
        for corner in 0..3 {
            let start = 1 + corner * 6;
            triangle[corner] = point(&record[start..start + 3])?;
            directions[corner] = normal(&record[start + 3..start + 6])?;
        }
        // Rounded native tessellation can collapse a face. Its area is exactly
        // zero; removing it changes neither the boundary nor its analytic normals.
        if crate::triangulation::triangle_cross(triangle[0], triangle[1], triangle[2])
            .iter()
            .all(|&value| value == 0.0)
        {
            discarded += 1;
            continue;
        }
        let mut face = [0; 3];
        for ((position, direction), index) in triangle.into_iter().zip(directions).zip(&mut face) {
            let key = (
                position.map(f32::to_bits),
                direction.map(f32::to_bits),
                atom,
            );
            *index = match vertices.entry(key) {
                std::collections::hash_map::Entry::Occupied(entry) => *entry.get(),
                std::collections::hash_map::Entry::Vacant(entry) => {
                    if positions.len() >= limits.vertices {
                        return Err(SurfaceOutputError::Limit("vertices"));
                    }
                    let row = u32::try_from(positions.len())
                        .map_err(|_| SurfaceOutputError::Limit("vertices"))?;
                    positions.push(position);
                    normals.push(direction);
                    atom_rows.push(atom);
                    entry.insert(row);
                    row
                }
            };
        }
        faces.push(SurfaceFace(face));
    }
    finish(
        positions,
        normals,
        atom_rows,
        faces,
        discarded,
        cancellation,
    )
}

fn finish(
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    atom_rows: Vec<u32>,
    faces: Vec<SurfaceFace>,
    discarded_degenerate_faces: u64,
    cancellation: &CancellationToken,
) -> Result<MolecularSurfaceMesh, SurfaceOutputError> {
    cancellation.check()?;
    if positions.is_empty() || faces.is_empty() {
        return Err(SurfaceOutputError::Invalid("surface contains no triangles"));
    }
    if faces.iter().any(|face| {
        face.0.iter().any(|&row| match usize::try_from(row) {
            Ok(row) => row >= positions.len(),
            Err(_) => true,
        })
    }) {
        return Err(SurfaceOutputError::Invalid("face addresses no vertex"));
    }
    let geometry =
        IndexedSurfaceMesh::from_analytic_normals(positions, faces, normals, cancellation)?;
    if geometry.report.degenerate_faces != 0 {
        return Err(SurfaceOutputError::Invalid(
            "surface contains degenerate faces",
        ));
    }
    Ok(MolecularSurfaceMesh {
        geometry,
        atom_rows: atom_rows.into_boxed_slice(),
        discarded_degenerate_faces,
    })
}

fn index(value: f64) -> Result<u32, SurfaceOutputError> {
    if value.fract() != 0.0 {
        return Err(SurfaceOutputError::Invalid("index is not an integer"));
    }
    value
        .to_u32()
        .ok_or(SurfaceOutputError::Invalid("index exceeds source rows"))
}

fn atom_row(rows: &[u32], index: u32) -> Result<u32, SurfaceOutputError> {
    rows.get(usize::try_from(index).map_err(|_| SurfaceOutputError::Limit("atom rows"))?)
        .copied()
        .ok_or(SurfaceOutputError::Invalid(
            "vertex addresses no input atom",
        ))
}

fn point(values: &[f64]) -> Result<[f32; 3], SurfaceOutputError> {
    let mut output = [0.0; 3];
    for (out, &value) in output.iter_mut().zip(values) {
        *out =
            value
                .to_f32()
                .filter(|value| value.is_finite())
                .ok_or(SurfaceOutputError::Invalid(
                    "coordinate exceeds finite precision",
                ))?;
    }
    Ok(output)
}

fn normal(values: &[f64]) -> Result<[f32; 3], SurfaceOutputError> {
    let length = values[0].hypot(values[1]).hypot(values[2]);
    if !length.is_finite() || length == 0.0 {
        return Err(SurfaceOutputError::Invalid(
            "vertex normal is not a finite direction",
        ));
    }
    point(&[values[0] / length, values[1] / length, values[2] / length])
}

struct Records<'a, R> {
    input: R,
    line: String,
    bytes: u64,
    maximum: u64,
    cancellation: &'a CancellationToken,
}

impl<'a, R: BufRead> Records<'a, R> {
    fn new(input: R, maximum: u64, cancellation: &'a CancellationToken) -> Self {
        Self {
            input,
            line: String::with_capacity(1025),
            bytes: 0,
            maximum,
            cancellation,
        }
    }

    fn next<const N: usize>(&mut self) -> Result<Option<[f64; N]>, SurfaceOutputError> {
        let mut record = [0.0; N];
        let mut fields = 0;
        loop {
            self.cancellation.check()?;
            self.line.clear();
            let read = (&mut self.input).take(1025).read_line(&mut self.line)?;
            if read == 0 {
                return if fields == 0 {
                    Ok(None)
                } else {
                    Err(SurfaceOutputError::Invalid("truncated record"))
                };
            }
            if read > 1024 {
                return Err(SurfaceOutputError::Limit("line length"));
            }
            self.bytes = self.bytes.saturating_add(
                u64::try_from(read).map_err(|_| SurfaceOutputError::Limit("bytes"))?,
            );
            if self.bytes > self.maximum {
                return Err(SurfaceOutputError::Limit("bytes"));
            }
            for token in self.line.split_ascii_whitespace() {
                let Some(out) = record.get_mut(fields) else {
                    return Err(SurfaceOutputError::Invalid("extra record fields"));
                };
                *out = token
                    .parse::<f64>()
                    .map_err(|_| SurfaceOutputError::Invalid("invalid number"))?;
                if !out.is_finite() {
                    return Err(SurfaceOutputError::Invalid("nonfinite number"));
                }
                fields += 1;
            }
            if fields == N {
                return Ok(Some(record));
            }
        }
    }
}

#[cfg(test)]
#[path = "solver_output_tests.rs"]
mod tests;
