//! Typed structure-local columns attached to atom rows.

use crate::column::{Presence, ValidityMask};
use crate::limits::CapacityError;
use crate::symbol::SymbolId;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Conventional custom column name for atom segment identifiers.
pub const SEGMENT_ID_ANNOTATION: &str = "segid";
/// Conventional custom column name for partial atomic charge.
pub const PARTIAL_CHARGE_ANNOTATION: &str = "charge";
/// Conventional custom column name for per-atom PQR radii.
pub const ATOM_RADIUS_ANNOTATION: &str = "radius";
/// Conventional custom column name for `AutoDock` atom types.
pub const AUTODOCK_TYPE_ANNOTATION: &str = "autodock_type";
/// Component role assigned from a chemical component definition.
pub const COMPONENT_KIND_ANNOTATION: &str = "component_kind";
/// Explicit polymer atom role assigned from a CCD-aware or caller-supplied profile.
pub const POLYMER_ATOM_ROLE_ANNOTATION: &str = "polymer_atom_role";
/// Aromatic-atom flag assigned from a chemical component definition.
pub const AROMATIC_ATOM_ANNOTATION: &str = "aromatic";
/// Formal charge assigned from a chemical component definition.
pub const FORMAL_CHARGE_ANNOTATION: &str = "formalcharge";
/// Hydrogen-bond donor role perceived from CCD connectivity and protonation.
pub const HBOND_DONOR_ANNOTATION: &str = "hbond_donor";
/// Hydrogen-bond acceptor role perceived from CCD connectivity and protonation.
pub const HBOND_ACCEPTOR_ANNOTATION: &str = "hbond_acceptor";
/// Modelled stereochemical configuration anchored to component coordinates.
pub const STEREO_CONFIGURATION_ANNOTATION: &str = "chirality";
/// Conventional custom column name for per-atom predicted confidence.
pub const PLDDT_ANNOTATION: &str = "plddt";
/// Conventional custom column name for per-atom predicted aligned error.
pub const PAE_ANNOTATION: &str = "pae";

/// One immutable typed annotation column with three-state validity.
#[derive(Clone, Debug)]
pub struct AnnotationColumn<T> {
    values: Arc<Vec<T>>,
    validity: ValidityMask,
}

impl<T> AnnotationColumn<T> {
    /// Builds an all-present column.
    ///
    /// # Errors
    ///
    /// Returns [`CapacityError`] when the column exceeds the supported row count.
    pub fn from_values(values: Vec<T>) -> Result<Self, CapacityError> {
        let len =
            u32::try_from(values.len()).map_err(|_| CapacityError::new("annotation column"))?;
        let validity = ValidityMask::all_present(len);
        Ok(Self {
            values: Arc::new(values),
            validity,
        })
    }

    /// Number of atom rows covered by the column.
    #[must_use]
    pub fn len(&self) -> u32 {
        self.validity.len()
    }

    /// Returns true when the column covers no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.validity.is_empty()
    }

    /// Validity state at one atom row.
    #[must_use]
    pub fn presence(&self, atom: u32) -> Presence {
        self.validity.get(atom)
    }

    /// Dense values backing this column.
    #[must_use]
    pub fn values(&self) -> &[T] {
        &self.values
    }

    pub(crate) fn filter(&self, keep: impl Fn(u32) -> bool) -> Result<Self, CapacityError>
    where
        T: Clone,
    {
        let entries = self
            .values
            .iter()
            .enumerate()
            .filter_map(|(position, value)| {
                let position = u32::try_from(position).ok()?;
                keep(position).then(|| (value.clone(), self.validity.get(position)))
            });
        Self::from_entries(entries)
    }
}

impl<T> AnnotationColumn<T> {
    /// Builds a column from dense values paired with their validity state.
    ///
    /// # Errors
    ///
    /// Returns [`CapacityError`] when the column exceeds the supported row count.
    pub fn from_entries(
        entries: impl IntoIterator<Item = (T, Presence)>,
    ) -> Result<Self, CapacityError> {
        let (values, presences): (Vec<_>, Vec<_>) = entries.into_iter().unzip();
        let validity = ValidityMask::try_from_iter(presences)
            .ok_or_else(|| CapacityError::new("annotation column"))?;
        Ok(Self {
            values: Arc::new(values),
            validity,
        })
    }
}

impl<T: Copy> AnnotationColumn<T> {
    /// Value and validity at one row, or `None` outside the column.
    #[must_use]
    pub fn get(&self, atom: u32) -> Option<(T, Presence)> {
        self.values
            .get(atom as usize)
            .copied()
            .map(|value| (value, self.validity.get(atom)))
    }
}

/// The supported physical types for a custom atom annotation.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum AtomAnnotation {
    /// Boolean values.
    Boolean(AnnotationColumn<bool>),
    /// Signed 64-bit integers.
    Integer(AnnotationColumn<i64>),
    /// IEEE-754 64-bit real values.
    Real(AnnotationColumn<f64>),
    /// Dictionary-backed categorical or textual values.
    Symbol(AnnotationColumn<SymbolId>),
}

impl AtomAnnotation {
    /// Number of atom rows covered by this column.
    #[must_use]
    pub fn len(&self) -> u32 {
        match self {
            Self::Boolean(column) => column.len(),
            Self::Integer(column) => column.len(),
            Self::Real(column) => column.len(),
            Self::Symbol(column) => column.len(),
        }
    }

    /// Returns true when the column covers no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub(crate) fn filter(&self, keep: impl Fn(u32) -> bool + Copy) -> Result<Self, CapacityError> {
        Ok(match self {
            Self::Boolean(column) => Self::Boolean(column.filter(keep)?),
            Self::Integer(column) => Self::Integer(column.filter(keep)?),
            Self::Real(column) => Self::Real(column.filter(keep)?),
            Self::Symbol(column) => Self::Symbol(column.filter(keep)?),
        })
    }
}

/// Custom per-atom columns, ordered by name for deterministic iteration.
#[derive(Clone, Debug, Default)]
pub struct AtomAnnotations {
    columns: BTreeMap<Box<str>, AtomAnnotation>,
}

impl AtomAnnotations {
    /// Number of custom columns.
    #[must_use]
    pub fn len(&self) -> usize {
        self.columns.len()
    }

    /// Returns true when no custom column is present.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }

    /// Looks up a column by its exact name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&AtomAnnotation> {
        self.columns.get(name)
    }

    /// Iterates columns in lexicographic name order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &AtomAnnotation)> {
        self.columns
            .iter()
            .map(|(name, column)| (name.as_ref(), column))
    }

    /// Inserts or replaces a column, returning the previous one.
    pub fn insert(
        &mut self,
        name: impl Into<Box<str>>,
        column: AtomAnnotation,
    ) -> Option<AtomAnnotation> {
        self.columns.insert(name.into(), column)
    }

    /// Removes a named column.
    pub fn remove(&mut self, name: &str) -> Option<AtomAnnotation> {
        self.columns.remove(name)
    }

    pub(crate) fn filter(&self, keep: impl Fn(u32) -> bool + Copy) -> Result<Self, CapacityError> {
        let columns = self
            .columns
            .iter()
            .map(|(name, column)| Ok((name.clone(), column.filter(keep)?)))
            .collect::<Result<_, CapacityError>>()?;
        Ok(Self { columns })
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
