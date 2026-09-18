//! Architecture-width non-periodic comparison kernel.

use super::brute::{allowed_pair, euclidean_distance_squared, finite, valid_position};
#[cfg(not(target_arch = "x86_64"))]
use wide::f32x4;
#[cfg(target_arch = "x86_64")]
use wide::f32x8;

/// Compares one valid left atom against all non-periodic right atoms.
///
/// Targets use eight lanes on x86-64 and four lanes on `AArch64`. Scalar tails
/// preserve identical emission semantics without temporary allocation.
#[cfg(not(target_arch = "x86_64"))]
pub(super) fn append_pairs<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left_atom: u32,
    left_position: [f32; 3],
    right: &[u32],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    let (batches, tail) = right.as_chunks::<4>();
    for atoms in batches {
        compare_four::<UNIQUE>(
            positions,
            left_atom,
            left_position,
            atoms,
            cutoff_squared,
            emit,
        );
    }
    append_tail::<UNIQUE>(
        positions,
        left_atom,
        left_position,
        tail,
        cutoff_squared,
        emit,
    );
}

#[cfg(target_arch = "x86_64")]
pub(super) fn append_pairs<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left_atom: u32,
    left_position: [f32; 3],
    right: &[u32],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    let (batches, tail) = right.as_chunks::<8>();
    for atoms in batches {
        compare_eight::<UNIQUE>(
            positions,
            left_atom,
            left_position,
            atoms,
            cutoff_squared,
            emit,
        );
    }
    append_tail::<UNIQUE>(
        positions,
        left_atom,
        left_position,
        tail,
        cutoff_squared,
        emit,
    );
}

#[cfg(not(target_arch = "x86_64"))]
fn compare_four<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left_atom: u32,
    left: [f32; 3],
    atoms: &[u32],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    let Ok(atom_ids) = <&[u32; 4]>::try_from(atoms) else {
        return;
    };
    let right = atom_ids.map(|atom| masked_position::<UNIQUE>(positions, left_atom, atom));
    let dx = f32x4::from(right.map(|position| position[0])) - f32x4::splat(left[0]);
    let dy = f32x4::from(right.map(|position| position[1])) - f32x4::splat(left[1]);
    let dz = f32x4::from(right.map(|position| position[2])) - f32x4::splat(left[2]);
    emit_lanes::<UNIQUE, 4>(
        left_atom,
        atom_ids,
        (dx * dx + dy * dy + dz * dz).to_array(),
        cutoff_squared,
        emit,
    );
}

#[cfg(target_arch = "x86_64")]
fn compare_eight<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left_atom: u32,
    left: [f32; 3],
    atoms: &[u32],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    let Ok(atom_ids) = <&[u32; 8]>::try_from(atoms) else {
        return;
    };
    let right = atom_ids.map(|atom| masked_position::<UNIQUE>(positions, left_atom, atom));
    let dx = f32x8::from(right.map(|position| position[0])) - f32x8::splat(left[0]);
    let dy = f32x8::from(right.map(|position| position[1])) - f32x8::splat(left[1]);
    let dz = f32x8::from(right.map(|position| position[2])) - f32x8::splat(left[2]);
    emit_lanes::<UNIQUE, 8>(
        left_atom,
        atom_ids,
        (dx * dx + dy * dy + dz * dz).to_array(),
        cutoff_squared,
        emit,
    );
}

fn masked_position<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left_atom: u32,
    atom: u32,
) -> [f32; 3] {
    const MASKED: [f32; 3] = [f32::NAN; 3];
    match usize::try_from(atom)
        .ok()
        .and_then(|index| positions.get(index))
        .copied()
        .filter(|position| finite(*position) && allowed_pair::<UNIQUE>(left_atom, atom))
    {
        Some(position) => position,
        None => MASKED,
    }
}

fn emit_lanes<const UNIQUE: bool, const LANES: usize>(
    left_atom: u32,
    atoms: &[u32; LANES],
    squared: [f32; LANES],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    for lane in 0..LANES {
        append_if_within::<UNIQUE>(left_atom, atoms[lane], squared[lane], cutoff_squared, emit);
    }
}

fn append_tail<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left_atom: u32,
    left: [f32; 3],
    atoms: &[u32],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    for &right_atom in atoms {
        if !allowed_pair::<UNIQUE>(left_atom, right_atom) {
            continue;
        }
        let Some(right) = valid_position(positions, right_atom) else {
            continue;
        };
        let squared = euclidean_distance_squared(left, right);
        append_if_within::<UNIQUE>(left_atom, right_atom, squared, cutoff_squared, emit);
    }
}

#[inline]
fn append_if_within<const UNIQUE: bool>(
    left_atom: u32,
    right_atom: u32,
    squared: f32,
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    if allowed_pair::<UNIQUE>(left_atom, right_atom) && squared <= cutoff_squared {
        emit(left_atom, right_atom, squared);
    }
}
