//! Resolving a file's description of its space group to exact operations.
//!
//! A file names its group three ways, in decreasing precision: a Hall symbol,
//! which fixes setting and origin; a Hermann–Mauguin symbol, which fixes the
//! setting except, for the rhombohedral types, which axes the cell is given on;
//! and an International Tables number, which fixes only the type. Each is used
//! for what it fixes and nothing more. The wrong setting of the right type is
//! the dangerous failure: it raises no error and places every symmetry mate in
//! the wrong place, so a symbol the catalogue cannot place is left unresolved
//! rather than replaced by a guess.

use super::{SpaceGroupSetting, catalogue_error, normalized, settings, symbol_key};
use crate::{Rational, SymmetryOperation, SymmetrySet};
use molframe_core::Diagnostic;
use molframe_core::structure::UnitCell;
use std::collections::BTreeSet;

/// Angular tolerance, in degrees, for reading a cell as rhombohedral.
const ANGLE_TOLERANCE: f64 = 0.01;

/// A Hall symbol resolved to its operations.
pub(super) enum HallResolution {
    /// The symbol, after any change of basis, is a catalogue setting.
    Setting(&'static SpaceGroupSetting),
    /// The change of basis leads outside the catalogue; the operations are
    /// still exact, but the setting has no catalogue number.
    Derived {
        international_number: u16,
        operations: Vec<SymmetryOperation>,
    },
}

/// Resolves a Hall symbol, including one with a change-of-basis suffix.
///
/// `P 61 2 (0 0 5)` and `P 61 2 (x,y,z+5/12)` name the same setting, and the
/// catalogue spells only the first. A suffix `V` turns every operation `S` of
/// the base symbol into `V S V⁻¹`; the result is matched, as a set of
/// operations modulo lattice translations, against every setting of the
/// type, so a recognised setting keeps its number.
///
/// Cost: one pass over the catalogue for the exact spelling, otherwise
/// `O(k · n log n)` for the `k` settings of the type with `n` operations each.
///
/// # Errors
///
/// Returns `E6018` when the base symbol is not in the catalogue or the
/// change of basis is malformed or not unimodular.
pub(super) fn resolve_hall(symbol: &str) -> Result<HallResolution, Diagnostic> {
    let requested = normalized(symbol);
    let catalogue = settings()?;
    if let Some(setting) = catalogue
        .iter()
        .find(|setting| normalized(&setting.hall_symbol) == requested)
    {
        return Ok(HallResolution::Setting(setting));
    }
    let (base, basis) = split_basis(&requested)?;
    let (international_number, base_operations) = base_operations(catalogue, base)?;
    let operations = base_operations
        .iter()
        .enumerate()
        .map(|(index, operation)| conjugate(&basis, operation, index))
        .collect::<Result<Vec<_>, _>>()?;
    let wanted = operation_key(&operations);
    let matched = catalogue.iter().find(|setting| {
        setting.international_number == international_number
            && operation_key(&setting.operations) == wanted
    });
    Ok(match matched {
        Some(setting) => HallResolution::Setting(setting),
        None => HallResolution::Derived {
            international_number,
            operations,
        },
    })
}

/// Looks up a Hermann–Mauguin symbol in the light of the cell it describes.
///
/// The rhombohedral types have two settings that share one symbol: hexagonal
/// axes (`a = b`, `γ = 120°`) and rhombohedral axes (`a = b = c`,
/// `α = β = γ`). A symbol that does not say which — no `H` lattice letter as
/// the PDB writes it and no `:H`/`:R` suffix — takes its axes from the cell;
/// without a cell it means the hexagonal axes, the setting a deposition uses.
///
/// # Errors
///
/// Returns `E6018` when no setting carries the symbol.
pub fn space_group_by_hermann_mauguin_in_cell(
    symbol: &str,
    cell: Option<&UnitCell>,
) -> Result<&'static SpaceGroupSetting, Diagnostic> {
    let setting = super::space_group_by_hermann_mauguin(symbol)?;
    if !is_rhombohedral_type(setting) || names_axes(symbol) {
        return Ok(setting);
    }
    let Some(cell) = cell else {
        return Ok(setting);
    };
    let wanted = if has_rhombohedral_axes(cell) {
        "R"
    } else {
        "H"
    };
    settings()?
        .iter()
        .find(|candidate| {
            candidate.international_number == setting.international_number
                && candidate.choice.as_ref() == wanted
        })
        .ok_or_else(catalogue_error)
}

/// Fills the operations of a set that names its group but lists none.
///
/// The Hall symbol decides when there is one; else the Hermann–Mauguin
/// symbol; else the number, in its standard setting (or, for a rhombohedral
/// type, the setting the cell's axes call for). A number that contradicts the
/// symbol is reported rather than overridden either way. A Hermann–Mauguin
/// symbol the catalogue does not know leaves a type with several settings
/// unresolved, since its number alone cannot say which setting the file is in.
pub(crate) fn resolve_metadata_operations(
    set: &mut SymmetrySet,
    cell: Option<&UnitCell>,
) -> Result<(), Diagnostic> {
    if let Some(hall) = set.hall.as_deref() {
        return match resolve_hall(hall)? {
            HallResolution::Setting(setting) => apply(set, setting),
            HallResolution::Derived {
                international_number,
                operations,
            } => {
                check_number(set, international_number)?;
                set.international_number = Some(international_number);
                set.operations = operations;
                Ok(())
            }
        };
    }
    if let Some(symbol) = set.hermann_mauguin.as_deref()
        && let Ok(setting) = space_group_by_hermann_mauguin_in_cell(symbol, cell)
    {
        return apply(set, setting);
    }
    let Some(number) = set.international_number else {
        return Ok(());
    };
    let candidates = super::space_group_settings(number)?;
    let setting = match candidates.as_slice() {
        [only] => *only,
        [first, ..] if set.hermann_mauguin.is_none() => {
            match cell.filter(|_| is_rhombohedral_type(first)) {
                Some(cell) => {
                    space_group_by_hermann_mauguin_in_cell(&first.international_full, Some(cell))?
                }
                None => *first,
            }
        }
        _ => return Ok(()),
    };
    apply(set, setting)
}

fn apply(set: &mut SymmetrySet, setting: &SpaceGroupSetting) -> Result<(), Diagnostic> {
    check_number(set, setting.international_number)?;
    set.hall_number = Some(setting.hall_number);
    set.international_number = Some(setting.international_number);
    if set.hermann_mauguin.is_none() {
        set.hermann_mauguin = Some(setting.international_short.clone());
    }
    if set.hall.is_none() {
        set.hall = Some(setting.hall_symbol.clone());
    }
    set.choice = (!setting.choice.is_empty()).then(|| setting.choice.clone());
    set.operations = setting.operations.to_vec();
    Ok(())
}

fn check_number(set: &SymmetrySet, resolved: u16) -> Result<(), Diagnostic> {
    match set.international_number {
        Some(declared) if declared != resolved => Err(catalogue_error()
            .with_context("declared_number", declared.to_string())
            .with_context("symbol_number", resolved.to_string())),
        _ => Ok(()),
    }
}

/// The operations of a bare Hall symbol, with its type number.
///
/// The catalogue spells some settings only with their own change of basis
/// (`P 61 2 (0 0 5)`), so a base symbol it lacks is recovered from such an
/// entry by undoing that entry's basis.
fn base_operations(
    catalogue: &'static [SpaceGroupSetting],
    base: &str,
) -> Result<(u16, Vec<SymmetryOperation>), Diagnostic> {
    if let Some(setting) = catalogue
        .iter()
        .find(|setting| normalized(&setting.hall_symbol) == base)
    {
        return Ok((setting.international_number, setting.operations.to_vec()));
    }
    for setting in catalogue {
        let symbol = normalized(&setting.hall_symbol);
        let Ok((entry_base, entry_basis)) = split_basis(&symbol) else {
            continue;
        };
        if entry_base != base {
            continue;
        }
        let undo = inverse(&entry_basis)?;
        let operations = setting
            .operations
            .iter()
            .enumerate()
            .map(|(index, operation)| conjugate(&undo, operation, index))
            .collect::<Result<Vec<_>, _>>()?;
        return Ok((setting.international_number, operations));
    }
    Err(catalogue_error())
}

/// The inverse affine map `(P⁻¹, −P⁻¹ p)` of a change of basis `(P, p)`.
fn inverse(basis: &SymmetryOperation) -> Result<SymmetryOperation, Diagnostic> {
    let rotation = unimodular_inverse(basis.rotation)?;
    let mut translation = [Rational::ZERO; 3];
    for (row, value) in translation.iter_mut().enumerate() {
        for (shift, factor) in basis.translation.iter().zip(rotation[row]) {
            *value = value.checked_sub(shift.checked_mul(factor)?)?;
        }
    }
    Ok(SymmetryOperation {
        id: "basis".into(),
        rotation,
        translation,
    })
}

/// Splits `base (basis)` and reads the basis as an operation: either an
/// algebraic triplet or three integer shifts in twelfths of the axes.
fn split_basis(symbol: &str) -> Result<(&str, SymmetryOperation), Diagnostic> {
    let inner = symbol.strip_suffix(')').ok_or_else(catalogue_error)?;
    let (base, basis) = inner.rsplit_once('(').ok_or_else(catalogue_error)?;
    let base = base.trim_end();
    if basis.contains(',') {
        let operation = SymmetryOperation::parse("basis", basis.trim())
            .map_err(|finding| catalogue_error().with_context("basis", finding.to_string()))?;
        return Ok((base, operation));
    }
    let shifts = basis
        .split_ascii_whitespace()
        .map(|shift| shift.parse::<i32>().map_err(|_| catalogue_error()))
        .collect::<Result<Vec<_>, _>>()?;
    let [x, y, z] = shifts.as_slice() else {
        return Err(catalogue_error());
    };
    Ok((
        base,
        SymmetryOperation {
            id: "basis".into(),
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: [
                Rational::new(*x, 12)?,
                Rational::new(*y, 12)?,
                Rational::new(*z, 12)?,
            ],
        },
    ))
}

/// `V S V⁻¹` for the basis `V = (P, p)` and operation `S = (W, w)`:
/// `W' = P W P⁻¹` and `w' = P w + p − W' p`, reduced into the unit cell.
fn conjugate(
    basis: &SymmetryOperation,
    operation: &SymmetryOperation,
    index: usize,
) -> Result<SymmetryOperation, Diagnostic> {
    let inverse = unimodular_inverse(basis.rotation)?;
    let rotation = multiply(multiply(basis.rotation, operation.rotation), inverse);
    let mut translation = [Rational::ZERO; 3];
    for (row, value) in translation.iter_mut().enumerate() {
        let mut sum = basis.translation[row];
        let terms = operation
            .translation
            .iter()
            .zip(basis.rotation[row])
            .zip(basis.translation.iter().zip(rotation[row]));
        for ((moved, into), (shift, factor)) in terms {
            sum = sum
                .checked_add(moved.checked_mul(into)?)?
                .checked_sub(shift.checked_mul(factor)?)?;
        }
        *value = reduced(sum)?;
    }
    Ok(SymmetryOperation {
        id: (index + 1).to_string().into_boxed_str(),
        rotation,
        translation,
    })
}

fn multiply(left: [[i32; 3]; 3], right: [[i32; 3]; 3]) -> [[i32; 3]; 3] {
    let mut product = [[0; 3]; 3];
    for (row, values) in product.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            *value = (0..3).map(|k| left[row][k] * right[k][column]).sum();
        }
    }
    product
}

/// The inverse of an integer matrix of determinant ±1, which is again integer.
fn unimodular_inverse(matrix: [[i32; 3]; 3]) -> Result<[[i32; 3]; 3], Diagnostic> {
    let m = matrix;
    let cofactor = |row: usize, column: usize| {
        let (r0, r1) = ((row + 1) % 3, (row + 2) % 3);
        let (c0, c1) = ((column + 1) % 3, (column + 2) % 3);
        m[r0][c0] * m[r1][c1] - m[r0][c1] * m[r1][c0]
    };
    let determinant: i32 = (0..3)
        .map(|column| m[0][column] * cofactor(0, column))
        .sum();
    if determinant.abs() != 1 {
        return Err(catalogue_error());
    }
    let mut inverse = [[0; 3]; 3];
    for (row, values) in inverse.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            *value = cofactor(column, row) * determinant;
        }
    }
    Ok(inverse)
}

/// A translation component moved into `[0, 1)`.
fn reduced(value: Rational) -> Result<Rational, Diagnostic> {
    let denominator = value.denominator();
    let modulus = i64::from(denominator);
    let numerator = i64::from(value.numerator()).rem_euclid(modulus);
    Rational::new(
        i32::try_from(numerator).map_err(|_| catalogue_error())?,
        denominator,
    )
}

/// An operation set as a comparable value, translations taken modulo the lattice.
fn operation_key(operations: &[SymmetryOperation]) -> BTreeSet<([[i32; 3]; 3], [Rational; 3])> {
    operations
        .iter()
        .filter_map(|operation| {
            let mut translation = [Rational::ZERO; 3];
            for (value, source) in translation.iter_mut().zip(operation.translation) {
                *value = reduced(source).ok()?;
            }
            Some((operation.rotation, translation))
        })
        .collect()
}

fn is_rhombohedral_type(setting: &SpaceGroupSetting) -> bool {
    matches!(setting.choice.as_ref(), "H" | "R")
}

/// Whether the symbol itself says which axes it is given on.
fn names_axes(symbol: &str) -> bool {
    let key = symbol_key(symbol);
    key.starts_with('H') || key.contains(':')
}

fn has_rhombohedral_axes(cell: &UnitCell) -> bool {
    let [alpha, beta, gamma] = cell.angles;
    (alpha - beta).abs() < ANGLE_TOLERANCE
        && (beta - gamma).abs() < ANGLE_TOLERANCE
        && (alpha - 90.0).abs() >= ANGLE_TOLERANCE
}
