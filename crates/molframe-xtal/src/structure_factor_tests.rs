use super::{Displacement, ScatteringSite, StructureFactorCalculator, StructureFactorError};
use crate::{CellTransform, GaussianFormFactor, SymmetryOperation};
use molframe_core::{Element, structure::UnitCell};
use std::f64::consts::PI;

fn element(symbol: &str) -> Element {
    match Element::from_symbol(symbol) {
        Some(element) => element,
        None => panic!("{symbol} is an element"),
    }
}

fn cubic(length: f64) -> CellTransform {
    match CellTransform::new(&UnitCell {
        lengths: [length; 3],
        angles: [90.0; 3],
    }) {
        Ok(cell) => cell,
        Err(error) => panic!("the cubic cell is valid: {error:?}"),
    }
}

fn operations(expressions: &[&str]) -> Vec<SymmetryOperation> {
    expressions
        .iter()
        .map(
            |expression| match SymmetryOperation::parse("op", expression) {
                Ok(operation) => operation,
                Err(error) => panic!("{expression} parses: {error:?}"),
            },
        )
        .collect()
}

fn site(symbol: &str, position: [f64; 3], displacement: Displacement) -> ScatteringSite {
    ScatteringSite {
        element: element(symbol),
        position,
        occupancy: 1.0,
        displacement,
    }
}

fn form(symbol: &str) -> GaussianFormFactor {
    match GaussianFormFactor::xray(element(symbol)) {
        Some(factor) => factor,
        None => panic!("{symbol} has a form factor"),
    }
}

fn calculate(
    operations: &[SymmetryOperation],
    sites: &[ScatteringSite],
    hkl: [i32; 3],
) -> super::Complex64 {
    let calculator = match StructureFactorCalculator::new(cubic(10.0), operations) {
        Ok(calculator) => calculator,
        Err(error) => panic!("{error}"),
    };
    match calculator.calculate(sites, hkl) {
        Ok(value) => value,
        Err(error) => panic!("{error}"),
    }
}

#[test]
fn one_atom_in_p1_scatters_with_its_form_factor_and_phase() {
    let p1 = operations(&["x,y,z"]);
    let atom = site("C", [0.25, 0.0, 0.0], Displacement::Isotropic(0.0));
    let amplitude = calculate(&p1, &[atom], [1, 0, 0]);
    let expected = form("C").value(0.25 / 100.0);
    assert!((amplitude.norm() - expected).abs() < 1e-12);
    assert!((amplitude.arg() - 2.0 * PI * 0.25).abs() < 1e-12);
}

#[test]
fn occupancy_scales_and_isotropic_b_damps_the_amplitude() {
    let p1 = operations(&["x,y,z"]);
    let mut atom = site("O", [0.0; 3], Displacement::Isotropic(20.0));
    atom.occupancy = 0.5;
    let stol2 = 0.25 * 3.0 / 100.0;
    let amplitude = calculate(&p1, &[atom], [1, 1, 1]);
    let expected = 0.5 * form("O").value(stol2) * (-20.0 * stol2).exp();
    assert!((amplitude.re - expected).abs() < 1e-12);
    assert!(amplitude.im.abs() < 1e-12);
}

#[test]
fn an_inversion_centre_makes_every_amplitude_real() {
    let p_one_bar = operations(&["x,y,z", "-x,-y,-z"]);
    let atom = site("N", [0.13, 0.31, 0.07], Displacement::Isotropic(5.0));
    for hkl in [[1, 2, 3], [0, 1, 4], [2, 0, 1]] {
        let amplitude = calculate(&p_one_bar, &[atom], hkl);
        assert!(amplitude.im.abs() < 1e-12, "{hkl:?}: {amplitude:?}");
    }
}

#[test]
fn an_isotropic_tensor_matches_the_equivalent_isotropic_b() {
    let p21 = operations(&["x,y,z", "-x,y+1/2,-z"]);
    let b = 12.0;
    let u = b / (8.0 * PI * PI);
    let isotropic = site("S", [0.2, 0.15, 0.35], Displacement::Isotropic(b));
    let tensor = site(
        "S",
        [0.2, 0.15, 0.35],
        Displacement::Anisotropic([u, u, u, 0.0, 0.0, 0.0]),
    );
    for hkl in [[1, 0, 0], [2, 3, 1], [0, 2, 5]] {
        let a = calculate(&p21, &[isotropic], hkl);
        let c = calculate(&p21, &[tensor], hkl);
        assert!((a.re - c.re).abs() < 1e-10 && (a.im - c.im).abs() < 1e-10);
    }
}

#[test]
fn a_missing_form_factor_and_an_empty_group_are_errors() {
    let p1 = operations(&["x,y,z"]);
    let calculator = match StructureFactorCalculator::new(cubic(10.0), &p1) {
        Ok(calculator) => calculator,
        Err(error) => panic!("{error}"),
    };
    let fermium = site("Fm", [0.0; 3], Displacement::Isotropic(0.0));
    assert_eq!(
        calculator.calculate(&[fermium], [1, 0, 0]),
        Err(StructureFactorError::MissingFormFactor(100))
    );
    assert_eq!(
        StructureFactorCalculator::new(cubic(10.0), &[]).err(),
        Some(StructureFactorError::NoOperations)
    );
}
