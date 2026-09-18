//! Charge-transfer iteration.

use super::parameters::{Coefficients, coefficients};
use super::perceive::component_inputs;
use super::{PeoeAtom, PeoeBond, PeoeError, PeoeOptions};
use crate::Component;

/// Charges every atom of a CCD component in component atom order.
///
/// # Errors
///
/// Returns [`PeoeError`] when component topology is inconsistent, an atomic
/// environment has no profile parameters, or options are invalid.
pub fn component_peoe_charges(
    component: &Component,
    options: PeoeOptions,
) -> Result<Vec<f64>, PeoeError> {
    let (atoms, bonds) = component_inputs(component)?;
    peoe_charges(&atoms, &bonds, options)
}

/// Runs PEOE on already perceived atom types without interpreting names.
///
/// # Errors
///
/// Returns [`PeoeError`] when a bond endpoint is invalid or options are invalid.
pub fn peoe_charges(
    atoms: &[PeoeAtom],
    bonds: &[PeoeBond],
    options: PeoeOptions,
) -> Result<Vec<f64>, PeoeError> {
    validate_options(options)?;
    validate_bonds(atoms.len(), bonds)?;
    let parameters: Vec<Coefficients> = atoms
        .iter()
        .map(|atom| coefficients(options.profile, atom.atom_type))
        .collect();
    let ionisation: Vec<f64> = atoms
        .iter()
        .zip(&parameters)
        .map(|(atom, values)| values.ionisation_electronegativity(atom.atom_type))
        .collect();
    let mut charges: Vec<f64> = atoms.iter().map(|atom| atom.formal_charge).collect();
    let mut electronegativities = vec![0.0; atoms.len()];
    let mut delta = vec![0.0; atoms.len()];
    let mut damping = options.initial_damping;
    for _ in 0..options.passes {
        update_electronegativities(&charges, &parameters, &mut electronegativities);
        delta.fill(0.0);
        transfer_bonds(
            bonds,
            &electronegativities,
            &ionisation,
            damping,
            options.minimum_electronegativity_difference,
            &mut delta,
        );
        for (charge, change) in charges.iter_mut().zip(&delta) {
            *charge += *change;
        }
        damping *= options.damping_factor;
    }
    Ok(charges)
}

fn update_electronegativities(charges: &[f64], parameters: &[Coefficients], output: &mut [f64]) {
    for ((charge, values), result) in charges.iter().zip(parameters).zip(output) {
        *result = values.a + *charge * (values.b + values.c * *charge);
    }
}

fn transfer_bonds(
    bonds: &[PeoeBond],
    electronegativities: &[f64],
    ionisation: &[f64],
    damping: f64,
    tolerance: f64,
    delta: &mut [f64],
) {
    for bond in bonds {
        let difference = electronegativities[bond.atom_b] - electronegativities[bond.atom_a];
        if difference.abs() <= tolerance {
            continue;
        }
        let denominator = if difference < 0.0 {
            ionisation[bond.atom_b]
        } else {
            ionisation[bond.atom_a]
        };
        let transfer = damping * difference / denominator;
        delta[bond.atom_a] += transfer;
        delta[bond.atom_b] -= transfer;
    }
}

fn validate_options(options: PeoeOptions) -> Result<(), PeoeError> {
    let valid = options.passes > 0
        && options.initial_damping.is_finite()
        && options.initial_damping > 0.0
        && options.damping_factor.is_finite()
        && options.damping_factor > 0.0
        && options.damping_factor <= 1.0
        && options.minimum_electronegativity_difference.is_finite()
        && options.minimum_electronegativity_difference >= 0.0;
    if valid {
        Ok(())
    } else {
        Err(PeoeError::InvalidOptions)
    }
}

fn validate_bonds(atom_count: usize, bonds: &[PeoeBond]) -> Result<(), PeoeError> {
    for bond in bonds {
        if bond.atom_a == bond.atom_b {
            return Err(PeoeError::SelfBond { atom: bond.atom_a });
        }
        if bond.atom_a >= atom_count {
            return Err(PeoeError::BondIndexOutsideAtomArray {
                atom: bond.atom_a,
                atom_count,
            });
        }
        if bond.atom_b >= atom_count {
            return Err(PeoeError::BondIndexOutsideAtomArray {
                atom: bond.atom_b,
                atom_count,
            });
        }
    }
    Ok(())
}
