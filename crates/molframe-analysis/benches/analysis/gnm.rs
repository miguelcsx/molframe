//! Dense reference and scalable fixture for the GNM benchmarks.

use molframe_analysis::GnmOptions;
use molframe_core::{ExecutionContext, selection::AtomSelection};
use molframe_spatial::pairs_within;
use nalgebra::{DMatrix, SymmetricEigen};

pub(super) fn dense_gnm_reference(
    positions: &[[f32; 3]],
    selection: &AtomSelection,
    options: GnmOptions,
) -> Vec<f64> {
    let selected: Vec<u32> = selection.into_iter().collect();
    let contacts = match pairs_within(
        positions,
        selection,
        selection,
        options.contact_distance,
        options.backend,
        None,
        &ExecutionContext::default(),
    ) {
        Ok(contacts) => contacts,
        Err(error) => panic!("dense GNM benchmark contact search failed: {error}"),
    };
    let mut matrix = DMatrix::zeros(selected.len(), selected.len());
    for contact in contacts {
        let (Ok(left), Ok(right)) = (
            selected.binary_search(&contact.first),
            selected.binary_search(&contact.second),
        ) else {
            panic!("dense GNM benchmark contact escaped its selection");
        };
        matrix[(left, right)] = -1.0;
        matrix[(right, left)] = -1.0;
        matrix[(left, left)] += 1.0;
        matrix[(right, right)] += 1.0;
    }
    let mut values: Vec<f64> = SymmetricEigen::new(matrix)
        .eigenvalues
        .iter()
        .copied()
        .filter(|value: &f64| value.abs() > options.zero_mode_tolerance)
        .collect();
    values.sort_by(f64::total_cmp);
    values.truncate(options.mode_count);
    values
}

pub(super) fn grid(x_count: u16, y_count: u16, z_count: u16) -> Vec<[f32; 3]> {
    let capacity = usize::from(x_count) * usize::from(y_count) * usize::from(z_count);
    let mut positions = Vec::with_capacity(capacity);
    for z in 0..z_count {
        for y in 0..y_count {
            for x in 0..x_count {
                positions.push([f32::from(x), f32::from(y), f32::from(z)]);
            }
        }
    }
    positions
}
