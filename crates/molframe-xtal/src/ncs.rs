//! Non-crystallographic symmetry operators and lazy generated images.

use crate::AffineTransform;
use crate::category_transform::{malformed, read_affine};
use molframe_cif::{Category, Document};
use molframe_core::{AtomIndex, Code, Diagnostic, ModelIndex, Structure};
use std::collections::BTreeMap;

/// Stable extension key used to attach NCS metadata to a structure.
pub const NCS_EXTENSION: &str = "molframe.xtal.ncs.v1";

/// Whether an NCS operator relates deposited coordinates or generates copies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NcsCode {
    /// All related coordinates are already present in the data block.
    Given,
    /// The operator generates coordinates absent from the data block.
    Generate,
}

/// One recorded non-crystallographic symmetry operator.
#[derive(Clone, Debug, PartialEq)]
pub struct NcsOperator {
    /// Source operator identifier.
    pub id: Box<str>,
    /// Dictionary-controlled application semantics.
    pub code: NcsCode,
    /// Optional source description.
    pub details: Option<Box<str>>,
    /// Recorded Cartesian affine transform.
    pub transform: AffineTransform,
}

/// Deterministically ordered NCS operator set.
#[derive(Clone, Debug, Default)]
pub struct NcsSet {
    operators: BTreeMap<Box<str>, NcsOperator>,
}

impl NcsSet {
    /// One operator by its source identifier.
    #[must_use]
    pub fn get(&self, id: &str) -> Option<&NcsOperator> {
        self.operators.get(id)
    }

    /// Every operator in lexical identifier order.
    pub fn operators(&self) -> impl Iterator<Item = &NcsOperator> {
        self.operators.values()
    }

    /// Operators that must generate coordinates absent from the input.
    pub fn generators(&self) -> impl Iterator<Item = &NcsOperator> {
        self.operators
            .values()
            .filter(|operator| operator.code == NcsCode::Generate)
    }

    /// Number of recorded operators.
    #[must_use]
    pub fn len(&self) -> usize {
        self.operators.len()
    }

    /// Whether no NCS operator was recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.operators.is_empty()
    }
}

/// Lowers `_struct_ncs_oper` from the first data block.
///
/// # Errors
///
/// Returns diagnostics for missing identifiers, invalid controlled codes,
/// duplicate identifiers, or incomplete/non-finite matrices.
pub fn lower_ncs(document: &Document) -> Result<NcsSet, Vec<Diagnostic>> {
    let Some(category) = document
        .first_block()
        .and_then(|block| block.category("struct_ncs_oper"))
    else {
        return Ok(NcsSet::default());
    };
    lower_category(category)
}

fn lower_category(category: &Category) -> Result<NcsSet, Vec<Diagnostic>> {
    let mut operators = BTreeMap::new();
    let mut findings = Vec::new();
    for row in 0..category.row_count() {
        let id = required_identifier(category, "id", row);
        let code = required_code(category, row);
        let transform = read_affine(category, row, Code::E6014);
        let (id, code, transform) = match (id, code, transform) {
            (Ok(id), Ok(code), Ok(transform)) => (id, code, transform),
            (id, code, transform) => {
                findings.extend(id.err());
                findings.extend(code.err());
                findings.extend(transform.err());
                continue;
            }
        };
        let operator = NcsOperator {
            id: id.clone(),
            code,
            details: category
                .identifier("details", row)
                .map(|value| value.into_owned().into_boxed_str()),
            transform,
        };
        if operators.insert(id.clone(), operator).is_some() {
            findings.push(malformed(category, row, "id", Code::E6014).with_context("value", id));
        }
    }
    if findings.is_empty() {
        Ok(NcsSet { operators })
    } else {
        Err(findings)
    }
}

fn required_identifier(
    category: &Category,
    item: &str,
    row: usize,
) -> Result<Box<str>, Diagnostic> {
    category
        .identifier(item, row)
        .filter(|value| !value.trim().is_empty())
        .map(|value| value.into_owned().into_boxed_str())
        .ok_or_else(|| malformed(category, row, item, Code::E6014))
}

fn required_code(category: &Category, row: usize) -> Result<NcsCode, Diagnostic> {
    let value = required_identifier(category, "code", row)?;
    if value.eq_ignore_ascii_case("given") {
        Ok(NcsCode::Given)
    } else if value.eq_ignore_ascii_case("generate") {
        Ok(NcsCode::Generate)
    } else {
        Err(malformed(category, row, "code", Code::E6014).with_context("value", value))
    }
}

/// One atom image generated lazily by an NCS operator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NcsAtomInstance<'a> {
    /// Atom in the deposited source structure.
    pub source_atom: AtomIndex,
    /// Operator that generates the image.
    pub operator: &'a NcsOperator,
}

/// Lazy generated NCS images over a structure snapshot.
#[derive(Clone, Debug)]
pub struct NcsView<'a> {
    structure: &'a Structure,
    set: &'a NcsSet,
}

impl NcsView<'_> {
    /// Number of generated whole-structure copies.
    #[must_use]
    pub fn copy_count(&self) -> usize {
        self.set.generators().count()
    }

    /// Generated atom/operator pairs without copied coordinates.
    pub fn atoms(&self) -> impl Iterator<Item = NcsAtomInstance<'_>> {
        self.set.generators().flat_map(move |operator| {
            (0..self.structure.atom_count()).map(move |atom| NcsAtomInstance {
                source_atom: AtomIndex::new(atom),
                operator,
            })
        })
    }

    /// Generated positions for one model, preserving source atom identity.
    pub fn positions(
        &self,
        model: ModelIndex,
    ) -> impl Iterator<Item = (NcsAtomInstance<'_>, Option<[f32; 3]>)> {
        let positions = self.structure.model_positions(model);
        self.atoms().map(move |instance| {
            let position = positions
                .and_then(|values| values.get(instance.source_atom.as_usize()))
                .copied()
                .map(|position| instance.operator.transform.apply(position));
            (instance, position)
        })
    }
}

/// Access to NCS metadata attached during structure reading.
pub trait NcsExt {
    /// Attached NCS operator set, when present.
    fn ncs_set(&self) -> Option<&NcsSet>;

    /// Lazy view containing only operators marked `generate`.
    fn ncs_generated(&self) -> Option<NcsView<'_>>;
}

impl NcsExt for Structure {
    fn ncs_set(&self) -> Option<&NcsSet> {
        self.extensions().get(NCS_EXTENSION)
    }

    fn ncs_generated(&self) -> Option<NcsView<'_>> {
        self.ncs_set().map(|set| NcsView {
            structure: self,
            set,
        })
    }
}

#[cfg(test)]
#[path = "ncs_tests.rs"]
mod tests;
