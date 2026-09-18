//! Discrete differential geometry of ordered backbone traces.
//!
//! Frames use only local Cα positions. Missing points and degenerate segments
//! invalidate the affected frame rather than fabricating an axis. Curvature and
//! torsion are expressed per ångström and are invariant under rigid motion.

use crate::numeric::exact_count;
use crate::{
    EigenError, EigenOptions, cross, dihedral, displacement, dot, gyration_axes_with_options, norm,
    normalise,
};

/// Local discrete Frenet data at one backbone position.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BackboneFrame {
    /// Unit tangent along the trace.
    pub tangent: [f64; 3],
    /// Unit normal toward the local turn.
    pub normal: [f64; 3],
    /// Unit binormal completing the right-handed frame.
    pub binormal: [f64; 3],
    /// Turning angle divided by mean adjacent segment length, in Å⁻¹.
    pub curvature: f64,
    /// Signed binormal rotation divided by mean local segment length, in Å⁻¹.
    pub torsion: Option<f64>,
}

/// Aggregate geometry of a contiguous helix-like trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelixGeometry {
    /// Unit principal axis, deterministically oriented from first to last point.
    pub axis: [f64; 3],
    /// Mean axial displacement per residue, in ångström.
    pub rise: f64,
    /// Mean signed rotation around the axis per residue, in radians.
    pub twist: f64,
}

/// Computes one optional Frenet frame per input position.
///
/// Endpoints have no centered local frame. A missing point invalidates only
/// frames whose local stencil includes it.
#[must_use]
pub fn backbone_frames(positions: &[Option<[f32; 3]>]) -> Vec<Option<BackboneFrame>> {
    let mut frames = vec![None; positions.len()];
    if positions.len() < 3 {
        return frames;
    }
    for centre in 1..positions.len() - 1 {
        let (Some(previous), Some(current), Some(next)) = (
            positions[centre - 1],
            positions[centre],
            positions[centre + 1],
        ) else {
            continue;
        };
        frames[centre] = local_frame(previous, current, next, torsion_at(positions, centre));
    }
    frames
}

/// Estimates axis, rise and twist for a contiguous trace.
///
/// # Errors
///
/// Returns `None` for fewer than four points, missing coordinates or a trace
/// without a stable axis.
pub fn helix_geometry(positions: &[Option<[f32; 3]>]) -> Result<Option<HelixGeometry>, EigenError> {
    helix_geometry_with_options(positions, EigenOptions::standard())
}

/// Estimates helix geometry with explicit eigensolver convergence controls.
///
/// # Errors
///
/// Returns [`EigenError`] when the principal axis cannot be decomposed.
pub fn helix_geometry_with_options(
    positions: &[Option<[f32; 3]>],
    options: EigenOptions,
) -> Result<Option<HelixGeometry>, EigenError> {
    if positions.len() < 4 || positions.iter().any(Option::is_none) {
        return Ok(None);
    }
    let points: Vec<_> = positions.iter().filter_map(|point| *point).collect();
    let (Some(first), Some(last)) = (points.first().copied(), points.last().copied()) else {
        return Ok(None);
    };
    let direction = displacement(first, last);
    let Some(decomposition) = gyration_axes_with_options(&points, options)? else {
        return Ok(None);
    };
    let mut axis = decomposition.dominant();
    if dot(axis, direction) < 0.0 {
        axis = axis.map(|value| -value);
    }
    let mut rise = 0.0;
    let mut twist = 0.0;
    let mut twist_count = 0usize;
    for pair in points.windows(2) {
        rise += dot(displacement(pair[0], pair[1]), axis);
    }
    let centre = centroid_f64(&points).ok_or(EigenError::InputTooLarge)?;
    for pair in points.windows(2) {
        let (Some(left), Some(right)) =
            (radial(pair[0], centre, axis), radial(pair[1], centre, axis))
        else {
            return Ok(None);
        };
        twist += signed_angle(left, right, axis);
        twist_count += 1;
    }
    if twist_count == 0 {
        return Ok(None);
    }
    let steps = exact_count(points.len() - 1).ok_or(EigenError::InputTooLarge)?;
    let twist_denominator = exact_count(twist_count).ok_or(EigenError::InputTooLarge)?;
    Ok(Some(HelixGeometry {
        axis,
        rise: rise / steps,
        twist: twist / twist_denominator,
    }))
}

fn local_frame(
    previous: [f32; 3],
    current: [f32; 3],
    next: [f32; 3],
    torsion_angle: Option<f64>,
) -> Option<BackboneFrame> {
    let incoming = displacement(previous, current);
    let outgoing = displacement(current, next);
    let incoming_length = norm(incoming);
    let outgoing_length = norm(outgoing);
    let tangent = normalise(displacement(previous, next))?;
    let binormal = normalise(cross(incoming, outgoing))?;
    let normal = normalise(cross(binormal, tangent))?;
    let cosine = (dot(incoming, outgoing) / (incoming_length * outgoing_length)).clamp(-1.0, 1.0);
    let mean_length = f64::midpoint(incoming_length, outgoing_length);
    if mean_length <= 0.0 {
        return None;
    }
    Some(BackboneFrame {
        tangent,
        normal,
        binormal,
        curvature: cosine.acos() / mean_length,
        torsion: torsion_angle.map(|angle| angle / mean_length),
    })
}

fn torsion_at(positions: &[Option<[f32; 3]>], centre: usize) -> Option<f64> {
    if centre < 2 || centre + 1 >= positions.len() {
        return None;
    }
    dihedral(
        positions[centre - 2]?,
        positions[centre - 1]?,
        positions[centre]?,
        positions[centre + 1]?,
    )
}

fn centroid_f64(points: &[[f32; 3]]) -> Option<[f64; 3]> {
    let count = exact_count(points.len())?;
    let mut centre = [0.0; 3];
    for point in points {
        for axis in 0..3 {
            centre[axis] += f64::from(point[axis]) / count;
        }
    }
    Some(centre)
}

fn radial(point: [f32; 3], centre: [f64; 3], axis: [f64; 3]) -> Option<[f64; 3]> {
    let offset = [
        f64::from(point[0]) - centre[0],
        f64::from(point[1]) - centre[1],
        f64::from(point[2]) - centre[2],
    ];
    let axial = dot(offset, axis);
    normalise([
        offset[0] - axial * axis[0],
        offset[1] - axial * axis[1],
        offset[2] - axial * axis[2],
    ])
}

fn signed_angle(left: [f64; 3], right: [f64; 3], axis: [f64; 3]) -> f64 {
    dot(axis, cross(left, right)).atan2(dot(left, right))
}

#[cfg(test)]
#[path = "geometry_tests.rs"]
mod tests;
