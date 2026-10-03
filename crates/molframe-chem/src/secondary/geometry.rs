//! Peptide continuity, reconstructed donors and backbone conformation.

use super::{DsspBackbone as Backbone, DsspOptions};

pub(super) fn squared_distance(left: Option<[f32; 3]>, right: Option<[f32; 3]>) -> f32 {
    let (Some(left), Some(right)) = (left, right) else {
        return f32::INFINITY;
    };
    left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
}

pub(super) fn linked(previous: &Backbone, next: &Backbone) -> bool {
    previous.model == next.model
        && previous.chain == next.chain
        && squared_distance(previous.carbon, next.nitrogen) <= 2.5 * 2.5
}

pub(super) fn continuous(backbones: &[Backbone], first: usize, last: usize) -> bool {
    backbones
        .get(first..=last)
        .is_some_and(|span| span.windows(2).all(|p| linked(&p[0], &p[1])))
}

pub(super) fn amide_hydrogen(
    backbones: &[Backbone],
    residue: usize,
    length: f32,
) -> Option<[f32; 3]> {
    let backbone = backbones.get(residue)?;
    if backbone.proline {
        return None;
    }
    let nitrogen = backbone.nitrogen?;
    let Some(previous) = residue.checked_sub(1).and_then(|i| backbones.get(i)) else {
        return Some(nitrogen);
    };
    // DSSP reconstructs from the preceding carbonyl even after a peptide gap;
    // continuity constrains local patterns, not this electrostatic estimate.
    if previous.model != backbone.model {
        return Some(nitrogen);
    }
    let (Some(carbon), Some(oxygen)) = (previous.carbon, previous.oxygen) else {
        return Some(nitrogen);
    };
    let norm = squared_distance(Some(carbon), Some(oxygen)).sqrt();
    if norm <= f32::EPSILON {
        return Some(nitrogen);
    }
    Some(core::array::from_fn(|axis| {
        nitrogen[axis] + (carbon[axis] - oxygen[axis]) * length / norm
    }))
}

pub(super) fn bond_energy(
    backbones: &[Backbone],
    acceptor: usize,
    donor: usize,
    options: &DsspOptions,
) -> f64 {
    let (Some(carbonyl), Some(amide)) = (backbones.get(acceptor), backbones.get(donor)) else {
        return 0.0;
    };
    let (Some(c), Some(o), Some(n), Some(h)) = (
        carbonyl.carbon,
        carbonyl.oxygen,
        amide.nitrogen,
        amide_hydrogen(backbones, donor, options.amide_hydrogen_distance),
    ) else {
        return 0.0;
    };
    let distances = [(o, n), (c, h), (o, h), (c, n)]
        .map(|(a, b)| f64::from(squared_distance(Some(a), Some(b))).sqrt());
    if distances.iter().any(|d| *d < 0.5) {
        return -9.9;
    }
    let [on, ch, oh, cn] = distances;
    let energy = options.electrostatic_prefactor * (1.0 / on + 1.0 / ch - 1.0 / oh - 1.0 / cn);
    ((energy * 1000.0).round() / 1000.0).max(-9.9)
}

pub(super) fn bends(backbones: &[Backbone], residue: usize, threshold: f32) -> bool {
    let Some(first) = residue.checked_sub(2) else {
        return false;
    };
    let last = residue + 2;
    if !continuous(backbones, first, last) {
        return false;
    }
    let (Some(before), Some(centre), Some(after)) = (
        backbones[first].ca,
        backbones[residue].ca,
        backbones[last].ca,
    ) else {
        return false;
    };
    let a = subtract(centre, before);
    let b = subtract(after, centre);
    let norms = (dot(a, a) * dot(b, b)).sqrt();
    norms > 0.0 && dot(a, b) < norms * f64::from(threshold).to_radians().cos()
}

fn subtract(a: [f32; 3], b: [f32; 3]) -> [f64; 3] {
    core::array::from_fn(|i| f64::from(a[i]) - f64::from(b[i]))
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn torsion(a: [f32; 3], b: [f32; 3], c: [f32; 3], d: [f32; 3]) -> Option<f64> {
    let ab = subtract(b, a);
    let bc = subtract(c, b);
    let cd = subtract(d, c);
    let first = cross(ab, bc);
    let second = cross(bc, cd);
    if dot(first, first) <= f64::EPSILON || dot(second, second) <= f64::EPSILON {
        return None;
    }
    Some(
        (dot(cross(first, second), bc) / dot(bc, bc).sqrt())
            .atan2(dot(first, second))
            .to_degrees(),
    )
}

pub(super) fn pp_conformation(backbones: &[Backbone], residue: usize) -> bool {
    let Some(before) = residue.checked_sub(1) else {
        return false;
    };
    if !continuous(backbones, before, residue + 1) {
        return false;
    }
    let current = &backbones[residue];
    let (Some(c0), Some(n), Some(ca), Some(c), Some(n1)) = (
        backbones[before].carbon,
        current.nitrogen,
        current.ca,
        current.carbon,
        backbones[residue + 1].nitrogen,
    ) else {
        return false;
    };
    matches!((torsion(c0,n,ca,c),torsion(n,ca,c,n1)), (Some(phi),Some(psi)) if (-104.0..=-46.0).contains(&phi) && (116.0..=174.0).contains(&psi))
}
