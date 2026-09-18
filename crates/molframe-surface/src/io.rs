//! Deterministic interchange writers for indexed surface meshes.

use crate::IndexedSurfaceMesh;
use std::io::{BufWriter, Write as _};
use std::path::Path;

/// Writes a triangular mesh as Wavefront OBJ without replacing an existing file.
///
/// Vertices, normals and faces retain their native deterministic ordering.
///
/// # Errors
///
/// Returns the underlying filesystem error if the destination exists or cannot
/// be written.
pub fn write_obj(path: impl AsRef<Path>, mesh: &IndexedSurfaceMesh) -> std::io::Result<()> {
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    let mut output = BufWriter::new(file);
    for vertex in &mesh.vertices {
        writeln!(output, "v {} {} {}", vertex[0], vertex[1], vertex[2])?;
    }
    for normal in &mesh.vertex_normals {
        writeln!(output, "vn {} {} {}", normal[0], normal[1], normal[2])?;
    }
    for face in &mesh.faces {
        let [first, second, third] = face.0.map(|index| index + 1);
        writeln!(
            output,
            "f {first}//{first} {second}//{second} {third}//{third}"
        )?;
    }
    output.flush()
}

#[cfg(test)]
#[path = "io_tests.rs"]
mod tests;
