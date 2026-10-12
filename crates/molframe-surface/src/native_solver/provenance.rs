//! Exact sphere validation for providers that mislabel isolated components.
//!
//! Only a complete analytical sphere can override a provider's atom index. Its
//! vertices must all lie on the original atom's sphere and its analytic normals
//! must face away from that center. Rounded textual output determines tolerance.
//! Cost is O(atoms + vertices) per component, with an explicit component-count ceiling.
use super::{NativeAtoms, NativeSurfaceError};
use crate::{MolecularSurfaceMesh, SurfaceOutputError};
use molframe_core::execution::ExecutionContext;

pub(super) fn validate_isolated_sphere(
    mesh: &mut MolecularSurfaceMesh,
    atoms: NativeAtoms<'_>,
    context: &ExecutionContext,
) -> Result<(), NativeSurfaceError> {
    let Some(&first) = mesh.atom_rows.first() else {
        return Ok(());
    };
    for chunk in mesh.atom_rows.chunks(1024) {
        context.cancellation().check()?;
        if chunk.iter().any(|&row| row != first) {
            return Ok(());
        }
    }
    if fits(mesh, atoms, first, context)? {
        return Ok(());
    }
    let mut matching = None;
    for (index, &row) in atoms.atoms.rows.iter().enumerate() {
        if index % 1024 == 0 {
            context.cancellation().check()?;
        }
        if matches_vertex(
            mesh.geometry.vertices[0],
            mesh.geometry.vertex_normals[0],
            atoms,
            row,
        ) {
            if matching.is_some() {
                return Err(
                    SurfaceOutputError::Invalid("ambiguous isolated sphere identity").into(),
                );
            }
            matching = Some(row);
        }
    }
    let Some(row) = matching else {
        return Err(
            SurfaceOutputError::Invalid("provider sphere identity contradicts geometry").into(),
        );
    };
    if !fits(mesh, atoms, row, context)? {
        return Err(SurfaceOutputError::Invalid(
            "provider component is not an original analytical sphere",
        )
        .into());
    }
    for chunk in mesh.atom_rows.chunks_mut(1024) {
        context.cancellation().check()?;
        chunk.fill(row);
    }
    Ok(())
}

fn fits(
    mesh: &MolecularSurfaceMesh,
    atoms: NativeAtoms<'_>,
    row: u32,
    context: &ExecutionContext,
) -> Result<bool, NativeSurfaceError> {
    for (index, (&point, &normal)) in mesh
        .geometry
        .vertices
        .iter()
        .zip(&mesh.geometry.vertex_normals)
        .enumerate()
    {
        if index.is_multiple_of(1024) {
            context.cancellation().check()?;
        }
        if !matches_vertex(point, normal, atoms, row) {
            return Ok(false);
        }
    }
    Ok(true)
}

fn matches_vertex(point: [f32; 3], normal: [f32; 3], atoms: NativeAtoms<'_>, row: u32) -> bool {
    let Some(center) = atoms.position(row) else {
        return false;
    };
    let radius = f64::from(atoms.atoms.radii[row as usize]);
    let displacement =
        std::array::from_fn::<_, 3, _>(|axis| f64::from(point[axis]) - f64::from(center[axis]));
    let length = displacement
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt();
    let alignment = displacement
        .iter()
        .zip(normal)
        .map(|(value, normal)| value * f64::from(normal))
        .sum::<f64>();
    (length - radius).abs() <= 0.002 && alignment >= length * 0.99999
}
