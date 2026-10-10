//! Watson–Crick edge and reference-geometry checks for candidate base pairs.
//!
//! Complementary identity plus a hydrogen-bond count is not enough to call a
//! pair canonical: a Hoogsteen or sugar-edge contact between an adenine and a
//! uracil also satisfies both. Two further tests separate them. The supporting
//! bond must join atoms on the Watson–Crick edges of both bases, and the two
//! bases must sit in the planar, roughly 10.5 Å wide arrangement a canonical
//! pair has.
//!
//! The edge atoms use the standard nucleobase atom names shared by the
//! canonical bases and by the modified bases that keep them.

use super::CanonicalBase;
use molframe_core::structure::ResidueRef;

/// Geometric limits that a Watson–Crick pair must satisfy.
///
/// The defaults are deliberately generous, so that distorted but genuine pairs
/// are kept while edge-on, stacked and non-planar contacts are not.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WatsonCrickGeometry {
    /// Smallest accepted C1'–C1' distance in ångström.
    pub minimum_c1_distance: f64,
    /// Largest accepted C1'–C1' distance in ångström. A canonical pair is about
    /// 10.5 Å wide, so the default window is 10.5 ± 1.5 Å.
    pub maximum_c1_distance: f64,
    /// Largest accepted angle in degrees between the two base-plane normals,
    /// which bounds the combined buckle and propeller twist. Normals are taken
    /// from the C2, C4 and C6 ring atoms, so the sense of the normal is
    /// ignored.
    pub maximum_plane_angle_degrees: f64,
}

impl Default for WatsonCrickGeometry {
    fn default() -> Self {
        Self {
            minimum_c1_distance: 9.0,
            maximum_c1_distance: 12.0,
            maximum_plane_angle_degrees: 40.0,
        }
    }
}

/// Whether a hydrogen bond between atoms `left` (of base `left_base`) and
/// `right` (of base `right_base`) lies on the Watson–Crick edges of both bases.
///
/// A–T/U pairs through N6…O4 and N1…N3; G–C pairs through O6…N4, N1…N3 and
/// N2…O2.
pub(super) fn watson_crick_contact(
    left_base: CanonicalBase,
    left: &str,
    right_base: CanonicalBase,
    right: &str,
) -> bool {
    let (first, second) = (left_base, right_base);
    let allowed: &[(&str, &str)] = match (first, second) {
        (CanonicalBase::Adenine, CanonicalBase::Thymine | CanonicalBase::Uracil) => {
            &[("N6", "O4"), ("N1", "N3")]
        }
        (CanonicalBase::Guanine, CanonicalBase::Cytosine) => {
            &[("O6", "N4"), ("N1", "N3"), ("N2", "O2")]
        }
        (CanonicalBase::Thymine | CanonicalBase::Uracil, CanonicalBase::Adenine) => {
            &[("O4", "N6"), ("N3", "N1")]
        }
        (CanonicalBase::Cytosine, CanonicalBase::Guanine) => {
            &[("N4", "O6"), ("N3", "N1"), ("O2", "N2")]
        }
        _ => return false,
    };
    allowed.contains(&(left, right))
}

/// Whether two residues have the reference geometry of a Watson–Crick pair.
///
/// Returns false when a required atom (C1', C2, C4 or C6) is missing, since the
/// geometry then cannot be shown to hold.
pub(super) fn has_watson_crick_geometry(
    first: ResidueRef<'_>,
    second: ResidueRef<'_>,
    limits: WatsonCrickGeometry,
) -> bool {
    let (Some(c1_first), Some(c1_second)) = (point(first, "C1'"), point(second, "C1'")) else {
        return false;
    };
    let distance = norm(sub(c1_first, c1_second));
    if distance < limits.minimum_c1_distance || distance > limits.maximum_c1_distance {
        return false;
    }
    let (Some(normal_first), Some(normal_second)) = (normal(first), normal(second)) else {
        return false;
    };
    let cosine = dot(normal_first, normal_second).abs();
    cosine >= limits.maximum_plane_angle_degrees.to_radians().cos()
}

fn point(residue: ResidueRef<'_>, name: &str) -> Option<[f64; 3]> {
    let [x, y, z] = residue.atom(name)?.position()?;
    Some([f64::from(x), f64::from(y), f64::from(z)])
}

fn normal(residue: ResidueRef<'_>) -> Option<[f64; 3]> {
    let c2 = point(residue, "C2")?;
    let c4 = point(residue, "C4")?;
    let c6 = point(residue, "C6")?;
    let cross = cross(sub(c4, c2), sub(c6, c2));
    let length = norm(cross);
    (length > 1e-9).then(|| [cross[0] / length, cross[1] / length, cross[2] / length])
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
