//! Bounded generation of atom images relevant to one crystal cutoff.

use crate::{CellTransform, SymmetrySet};
use pdbiox_core::structure::AtomRef;
use pdbiox_core::{AtomIndex, Code, Diagnostic, ModelIndex, Structure};

#[derive(Clone, Copy, Debug)]
pub(crate) struct SourcePosition {
    pub(crate) atom: AtomIndex,
    pub(crate) cartesian: [f32; 3],
    fractional: [f64; 3],
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct CrystalImage {
    pub(crate) atom: AtomIndex,
    pub(crate) operation: usize,
    pub(crate) lattice: [i32; 3],
    pub(crate) cartesian: [f32; 3],
}

pub(crate) fn relevant_images(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
    limit: usize,
) -> Result<(Vec<SourcePosition>, Vec<CrystalImage>), Diagnostic> {
    validate(structure, symmetry, cutoff)?;
    let cell = structure
        .data()
        .cell
        .ok_or_else(|| Diagnostic::new(Code::E5004))?;
    let transform = CellTransform::new(&cell)?;
    let sources = valid_positions(structure, model, &transform)?;
    if sources.is_empty() {
        return Err(Diagnostic::new(Code::E5003));
    }
    let (minimum, maximum) = fractional_extent(&sources);
    let bounds = fractional_bounds(&transform, cutoff);
    let mut images = Vec::new();
    let mut candidates = 0usize;
    for (operation_index, operation) in symmetry.operations().iter().enumerate() {
        for source in &sources {
            let transformed = operation.apply_fractional(source.fractional);
            let ranges = image_ranges(transformed, minimum, maximum, bounds)?;
            for lattice_x in ranges[0].clone() {
                for lattice_y in ranges[1].clone() {
                    for lattice_z in ranges[2].clone() {
                        candidates = candidates
                            .checked_add(sources.len())
                            .ok_or_else(search_limit)?;
                        if candidates > limit {
                            return Err(search_limit().with_context("limit", limit.to_string()));
                        }
                        let lattice = [lattice_x, lattice_y, lattice_z];
                        let fractional = [
                            transformed[0] + f64::from(lattice_x),
                            transformed[1] + f64::from(lattice_y),
                            transformed[2] + f64::from(lattice_z),
                        ];
                        images.push(CrystalImage {
                            atom: source.atom,
                            operation: operation_index,
                            lattice,
                            cartesian: transform.to_cartesian(fractional).map(|value| value as f32),
                        });
                    }
                }
            }
        }
    }
    Ok((sources, images))
}

fn validate(structure: &Structure, symmetry: &SymmetrySet, cutoff: f64) -> Result<(), Diagnostic> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(Diagnostic::new(Code::E6016).with_context("cutoff", cutoff.to_string()));
    }
    if symmetry.operations().is_empty() {
        return Err(Diagnostic::new(Code::E6016).with_context("symmetry", "no explicit operators"));
    }
    if structure.data().cell.is_none() {
        return Err(Diagnostic::new(Code::E5004));
    }
    Ok(())
}

fn valid_positions(
    structure: &Structure,
    model: ModelIndex,
    transform: &CellTransform,
) -> Result<Vec<SourcePosition>, Diagnostic> {
    let frame = structure
        .model_positions(model)
        .ok_or_else(|| Diagnostic::new(Code::E6003).with_context("model", model.to_string()))?;
    let mut output = Vec::new();
    for atom in 0..structure.atom_count() {
        let index = AtomIndex::new(atom);
        if structure
            .data()
            .atom(index)
            .and_then(AtomRef::position)
            .is_none()
        {
            continue;
        }
        let cartesian = *frame.get(atom as usize).ok_or_else(invariant)?;
        if !cartesian.iter().all(|value| value.is_finite()) {
            continue;
        }
        output.push(SourcePosition {
            atom: index,
            cartesian,
            fractional: transform.to_fractional(cartesian.map(f64::from)),
        });
    }
    Ok(output)
}

fn fractional_extent(sources: &[SourcePosition]) -> ([f64; 3], [f64; 3]) {
    let mut minimum = [f64::INFINITY; 3];
    let mut maximum = [f64::NEG_INFINITY; 3];
    for source in sources {
        for axis in 0..3 {
            minimum[axis] = minimum[axis].min(source.fractional[axis]);
            maximum[axis] = maximum[axis].max(source.fractional[axis]);
        }
    }
    (minimum, maximum)
}

fn fractional_bounds(transform: &CellTransform, cutoff: f64) -> [f64; 3] {
    transform
        .inverse_matrix()
        .map(|row| cutoff * row.iter().map(|value| value * value).sum::<f64>().sqrt())
}

fn image_ranges(
    image: [f64; 3],
    minimum: [f64; 3],
    maximum: [f64; 3],
    bounds: [f64; 3],
) -> Result<[std::ops::RangeInclusive<i32>; 3], Diagnostic> {
    let mut ranges = [0..=0, 0..=0, 0..=0];
    for axis in 0..3 {
        let start = (minimum[axis] - bounds[axis] - image[axis]).ceil();
        let end = (maximum[axis] + bounds[axis] - image[axis]).floor();
        if start < f64::from(i32::MIN) || end > f64::from(i32::MAX) {
            return Err(search_limit());
        }
        ranges[axis] = start as i32..=end as i32;
    }
    Ok(ranges)
}

fn search_limit() -> Diagnostic {
    Diagnostic::new(Code::E6017)
}

fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}
