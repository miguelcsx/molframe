//! Exact crystallographic symmetry operations from CIF algebraic notation.

use pdbiox_cif::{Category, DataBlock, Document};
use pdbiox_core::{Code, Diagnostic, ModelIndex, Structure};
use std::collections::BTreeSet;

#[path = "symmetry_parse.rs"]
mod parse;

use parse::{determinant, gcd, parse_component};

/// Stable extension key for crystallographic space-group metadata.
pub const SYMMETRY_EXTENSION: &str = "pdbiox.xtal.symmetry.v1";

/// An exact reduced rational number.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rational {
    numerator: i32,
    denominator: u32,
}

impl Rational {
    /// Zero.
    pub const ZERO: Self = Self {
        numerator: 0,
        denominator: 1,
    };

    /// Creates and reduces an exact rational.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for a zero denominator or unrepresentable result.
    pub fn new(numerator: i32, denominator: u32) -> Result<Self, Diagnostic> {
        if denominator == 0 {
            return Err(symmetry_error("zero rational denominator"));
        }
        let divisor = gcd(numerator.unsigned_abs(), denominator);
        let reduced = i64::from(numerator) / i64::from(divisor);
        Ok(Self {
            numerator: i32::try_from(reduced).map_err(|_| symmetry_capacity())?,
            denominator: denominator / divisor,
        })
    }

    /// Signed numerator of the reduced fraction.
    #[must_use]
    pub const fn numerator(self) -> i32 {
        self.numerator
    }

    /// Positive denominator of the reduced fraction.
    #[must_use]
    pub const fn denominator(self) -> u32 {
        self.denominator
    }

    /// Floating representation for coordinate arithmetic.
    #[must_use]
    pub fn as_f64(self) -> f64 {
        f64::from(self.numerator) / f64::from(self.denominator)
    }

    pub(crate) fn checked_add(self, other: Self) -> Result<Self, Diagnostic> {
        let left = i64::from(self.numerator) * i64::from(other.denominator);
        let right = i64::from(other.numerator) * i64::from(self.denominator);
        let denominator = u64::from(self.denominator) * u64::from(other.denominator);
        Self::new(
            i32::try_from(left + right).map_err(|_| symmetry_capacity())?,
            u32::try_from(denominator).map_err(|_| symmetry_capacity())?,
        )
    }

    pub(crate) fn checked_sub(self, other: Self) -> Result<Self, Diagnostic> {
        self.checked_add(other.negated()?)
    }

    pub(crate) fn checked_mul(self, factor: i32) -> Result<Self, Diagnostic> {
        let numerator = i64::from(self.numerator) * i64::from(factor);
        Self::new(
            i32::try_from(numerator).map_err(|_| symmetry_capacity())?,
            self.denominator,
        )
    }

    pub(crate) fn integer(self) -> Option<i32> {
        (self.denominator == 1).then_some(self.numerator)
    }

    fn negated(self) -> Result<Self, Diagnostic> {
        Self::new(
            self.numerator.checked_neg().ok_or_else(symmetry_capacity)?,
            self.denominator,
        )
    }
}

/// One exact fractional-coordinate space-group operation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SymmetryOperation {
    /// Source operation identifier.
    pub id: Box<str>,
    /// Integer rotational matrix `W` in `X' = W X + w`.
    pub rotation: [[i32; 3]; 3],
    /// Exact fractional translation `w`.
    pub translation: [Rational; 3],
}

impl SymmetryOperation {
    /// Parses a CIF algebraic operation such as `-x,1/2+y,1/2-z`.
    ///
    /// # Errors
    ///
    /// Returns `E6015` for invalid syntax or a non-unimodular rotation.
    pub fn parse(id: impl Into<Box<str>>, expression: &str) -> Result<Self, Diagnostic> {
        let components = expression.split(',').collect::<Vec<_>>();
        if components.len() != 3 {
            return Err(symmetry_error(expression));
        }
        let mut rotation = [[0; 3]; 3];
        let mut translation = [Rational::ZERO; 3];
        for (row, component) in components.iter().enumerate() {
            parse_component(component, &mut rotation[row], &mut translation[row])?;
        }
        if determinant(rotation).unsigned_abs() != 1 {
            return Err(symmetry_error(expression));
        }
        Ok(Self {
            id: id.into(),
            rotation,
            translation,
        })
    }

    /// Applies the operation to one fractional coordinate.
    #[must_use]
    pub fn apply_fractional(&self, position: [f64; 3]) -> [f64; 3] {
        let mut output = [0.0; 3];
        for (row, value) in output.iter_mut().enumerate() {
            *value = (0..3)
                .map(|axis| f64::from(self.rotation[row][axis]) * position[axis])
                .sum::<f64>()
                + self.translation[row].as_f64();
        }
        output
    }

    /// Whether this is the identity representative with no lattice shift.
    #[must_use]
    pub fn is_identity(&self) -> bool {
        self.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            && self.translation == [Rational::ZERO; 3]
    }
}

/// Space-group identity and explicit coordinate representatives.
#[derive(Clone, Debug, Default)]
pub struct SymmetrySet {
    /// Catalogue serial number for the Hall setting, when resolved.
    pub hall_number: Option<u16>,
    /// International Tables space-group number.
    pub international_number: Option<u16>,
    /// Hermann–Mauguin symbol in the file's setting.
    pub hermann_mauguin: Option<Box<str>>,
    /// Hall symbol, which uniquely identifies setting and origin.
    pub hall: Option<Box<str>>,
    /// Declared crystal system.
    pub crystal_system: Option<Box<str>>,
    /// Unique-axis, origin or cell choice within the space-group type.
    pub choice: Option<Box<str>>,
    pub(crate) operations: Vec<SymmetryOperation>,
}

impl SymmetrySet {
    /// Explicit representatives in file order.
    #[must_use]
    pub fn operations(&self) -> &[SymmetryOperation] {
        &self.operations
    }

    /// Whether neither metadata nor explicit operations were recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.international_number.is_none()
            && self.hall_number.is_none()
            && self.hermann_mauguin.is_none()
            && self.hall.is_none()
            && self.crystal_system.is_none()
            && self.choice.is_none()
            && self.operations.is_empty()
    }
}

/// Lowers modern and legacy CIF symmetry categories from the first block.
///
/// # Errors
///
/// Returns diagnostics for invalid group numbers, duplicate operation IDs, or
/// malformed algebraic operations.
pub fn lower_symmetry(document: &Document) -> Result<SymmetrySet, Vec<Diagnostic>> {
    let Some(block) = document.first_block() else {
        return Ok(SymmetrySet::default());
    };
    let mut set = metadata(block)?;
    let category = block
        .category("space_group_symop")
        .map(|category| (category, "operation_xyz"))
        .or_else(|| {
            block
                .category("symmetry_equiv")
                .map(|category| (category, "pos_as_xyz"))
        });
    if let Some((category, item)) = category {
        set.operations = lower_operations(category, item)?;
    }
    if set.operations.is_empty() {
        crate::space_group::resolve_metadata_operations(&mut set)
            .map_err(|finding| vec![finding])?;
    }
    Ok(set)
}

/// Access to crystallographic symmetry metadata attached during reading.
pub trait SymmetryExt {
    /// Attached space-group metadata and explicit operations, when present.
    fn symmetry_set(&self) -> Option<&SymmetrySet>;

    /// Finds unique crystal contacts in the first coordinate model.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when symmetry, cell, coordinates or cutoff are invalid.
    fn collect_crystal_neighbors(
        &self,
        cutoff: f64,
        context: &pdbiox_core::ExecutionContext,
    ) -> Result<Vec<crate::CrystalNeighbor>, Diagnostic>;

    /// Finds unique crystal contacts in a selected coordinate model.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when symmetry, cell, model or cutoff are invalid.
    fn collect_crystal_neighbors_in_model(
        &self,
        model: ModelIndex,
        cutoff: f64,
        context: &pdbiox_core::ExecutionContext,
    ) -> Result<Vec<crate::CrystalNeighbor>, Diagnostic>;
}

impl SymmetryExt for Structure {
    fn symmetry_set(&self) -> Option<&SymmetrySet> {
        self.extensions().get(SYMMETRY_EXTENSION)
    }

    fn collect_crystal_neighbors(
        &self,
        cutoff: f64,
        context: &pdbiox_core::ExecutionContext,
    ) -> Result<Vec<crate::CrystalNeighbor>, Diagnostic> {
        self.collect_crystal_neighbors_in_model(ModelIndex::new(0), cutoff, context)
    }

    fn collect_crystal_neighbors_in_model(
        &self,
        model: ModelIndex,
        cutoff: f64,
        context: &pdbiox_core::ExecutionContext,
    ) -> Result<Vec<crate::CrystalNeighbor>, Diagnostic> {
        let symmetry = self
            .symmetry_set()
            .ok_or_else(|| Diagnostic::new(Code::E6016).with_context("symmetry", "absent"))?;
        crate::collect_crystal_neighbors(
            self,
            symmetry,
            model,
            cutoff,
            crate::CrystalNeighborOptions::default(),
            context,
        )
    }
}

fn metadata(block: &DataBlock) -> Result<SymmetrySet, Vec<Diagnostic>> {
    if let Some(category) = block.category("space_group") {
        return modern_metadata(category);
    }
    let Some(category) = block.category("symmetry") else {
        return Ok(SymmetrySet::default());
    };
    let number = group_number(category, "Int_Tables_number")?;
    Ok(SymmetrySet {
        hall_number: None,
        international_number: number,
        hermann_mauguin: text(category, "space_group_name_H-M"),
        hall: text(category, "space_group_name_Hall"),
        crystal_system: text(category, "cell_setting"),
        choice: None,
        operations: Vec::new(),
    })
}

fn modern_metadata(category: &Category) -> Result<SymmetrySet, Vec<Diagnostic>> {
    Ok(SymmetrySet {
        hall_number: None,
        international_number: group_number(category, "IT_number")?,
        hermann_mauguin: text(category, "name_H-M_alt"),
        hall: text(category, "name_Hall"),
        crystal_system: text(category, "crystal_system"),
        choice: None,
        operations: Vec::new(),
    })
}

fn group_number(category: &Category, item: &str) -> Result<Option<u16>, Vec<Diagnostic>> {
    let Some(value) = category.value(item, 0) else {
        return Ok(None);
    };
    if !value.is_recorded() {
        return Ok(None);
    }
    value
        .as_integer()
        .and_then(|number| u16::try_from(number).ok())
        .filter(|number| (1..=230).contains(number))
        .map(Some)
        .ok_or_else(|| vec![symmetry_error(item)])
}

fn lower_operations(
    category: &Category,
    item: &str,
) -> Result<Vec<SymmetryOperation>, Vec<Diagnostic>> {
    let mut operations = Vec::with_capacity(category.row_count());
    let mut identifiers = BTreeSet::new();
    let mut findings = Vec::new();
    for row in 0..category.row_count() {
        let Some(id) = category.identifier("id", row) else {
            findings.push(symmetry_row_error(category, "id", row));
            continue;
        };
        let id = id.into_owned().into_boxed_str();
        let Some(expression) = category.identifier(item, row) else {
            findings.push(symmetry_row_error(category, item, row));
            continue;
        };
        if !identifiers.insert(id.clone()) {
            findings.push(symmetry_row_error(category, "id", row).with_context("value", id));
            continue;
        }
        match SymmetryOperation::parse(id, &expression) {
            Ok(operation) => operations.push(operation),
            Err(finding) => findings.push(
                finding
                    .with_context("category", category.name())
                    .with_context("row", row.to_string()),
            ),
        }
    }
    if findings.is_empty() {
        Ok(operations)
    } else {
        Err(findings)
    }
}

fn text(category: &Category, item: &str) -> Option<Box<str>> {
    category
        .identifier(item, 0)
        .map(|value| value.into_owned().into_boxed_str())
}

fn symmetry_row_error(category: &Category, item: &str, row: usize) -> Diagnostic {
    symmetry_error(item)
        .with_context("category", category.name())
        .with_context("row", row.to_string())
}

pub(crate) fn symmetry_error(value: &str) -> Diagnostic {
    Diagnostic::new(Code::E6015).with_context("value", value)
}

pub(crate) fn symmetry_capacity() -> Diagnostic {
    Diagnostic::new(Code::E1901).with_context("limit", "symmetry expression arithmetic")
}

#[cfg(test)]
#[path = "symmetry_tests.rs"]
mod tests;
