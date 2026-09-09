//! Allocation-free typed projection shared by canonical format writers.
//!
//! Validation is linear in emitted atom rows. A subsequent visit is also
//! linear and allocates nothing, so a columnar writer can replay one column at
//! a time without retaining a row DOM.

use super::options::{CifWriteError, CifWriteOptions, valid_block_id};
use pdbiox_core::element::Element;
use pdbiox_core::index::ModelIndex;
use pdbiox_core::structure::{
    AtomRef, ChainRef, ResidueRef, SEQUENCE_REFERENCES_EXTENSION, SequenceReferences, Structure,
};

/// A typed canonical value with CIF's two distinct missing-value states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CanonicalValue<T> {
    /// A recorded value.
    Present(T),
    /// The item does not apply to this row.
    Inapplicable,
    /// The item applies but was not recorded.
    Unknown,
}

/// A validated semantic projection that can be consumed without text or a DOM.
#[derive(Clone, Copy, Debug)]
pub struct CanonicalProjection<'a> {
    structure: &'a Structure,
    options: &'a CifWriteOptions,
    block_id: &'a str,
}

/// Validates and borrows the canonical structure projection.
///
/// # Errors
///
/// Returns a canonical write refusal before a format writer allocates its
/// output when a required model, atom, bond endpoint, or identifier is absent.
pub fn canonical_projection<'a>(
    structure: &'a Structure,
    options: &'a CifWriteOptions,
) -> Result<CanonicalProjection<'a>, CifWriteError> {
    let Some(block_id) = options.block_id().or(structure.data().entry.id.as_deref()) else {
        return Err(CifWriteError::MissingBlockId);
    };
    if !valid_block_id(block_id) {
        return Err(CifWriteError::InvalidBlockId(block_id.to_owned()));
    }
    let projection = CanonicalProjection {
        structure,
        options,
        block_id,
    };
    projection.atom_site_row_count()?;
    projection.visit_atom_rows(|_| {})?;
    super::bonds::preflight(structure, options)?;
    Ok(projection)
}

impl<'a> CanonicalProjection<'a> {
    /// The validated data-block identifier.
    #[must_use]
    pub const fn block_id(self) -> &'a str {
        self.block_id
    }

    /// The immutable semantic structure behind this projection.
    #[must_use]
    pub const fn structure(self) -> &'a Structure {
        self.structure
    }

    /// The explicit canonical decisions used for this projection.
    #[must_use]
    pub const fn options(self) -> &'a CifWriteOptions {
        self.options
    }

    /// Connection type explicitly chosen for generated canonical bond rows.
    #[must_use]
    pub fn connection_type_id(self) -> Option<&'a str> {
        self.options.connection_type_id()
    }

    /// Database sequence references retained by the structure.
    #[must_use]
    pub fn sequence_references(self) -> Option<&'a SequenceReferences> {
        self.structure
            .extensions()
            .get::<SequenceReferences>(SEQUENCE_REFERENCES_EXTENSION)
    }

    /// Number of rows emitted by the atom-site category.
    ///
    /// # Errors
    ///
    /// Returns a range refusal if the product or ragged-model sum exceeds
    /// `usize`.
    pub fn atom_site_row_count(self) -> Result<usize, CifWriteError> {
        if let Some(models) = self.structure.ragged_models() {
            return models.iter().try_fold(0usize, |total, model| {
                let rows = usize::try_from(model.atom_count())
                    .map_err(|_| CifWriteError::AtomRowCountOverflow)?;
                total
                    .checked_add(rows)
                    .ok_or(CifWriteError::AtomRowCountOverflow)
            });
        }
        let atoms = usize::try_from(self.structure.atom_count())
            .map_err(|_| CifWriteError::AtomRowCountOverflow)?;
        atoms
            .checked_mul(self.structure.model_count())
            .ok_or(CifWriteError::AtomRowCountOverflow)
    }

    /// Visits canonical atom-site rows in deposition and model order.
    ///
    /// The row borrows a model snapshot and is valid only for the synchronous
    /// callback. Replaying the visit is allocation-free and supports encoders
    /// that materialise only one output column at a time.
    ///
    /// # Errors
    ///
    /// Returns a canonical write refusal if the immutable structure cannot
    /// supply a required row identity.
    pub fn visit_atom_rows(
        self,
        mut visit: impl for<'row> FnMut(CanonicalAtomRow<'row>),
    ) -> Result<(), CifWriteError> {
        for position in 0..self.structure.model_count() {
            let model = model_index(position)?;
            let Some((snapshot, local_model)) = self.structure.model_snapshot(model) else {
                return Err(CifWriteError::MissingModelNumber { model: model.get() });
            };
            let Some(model_number) = snapshot
                .model(local_model)
                .and_then(pdbiox_core::structure::ModelRef::number)
            else {
                return Err(CifWriteError::MissingModelNumber { model: model.get() });
            };
            visit_model(&snapshot, local_model, model_number, &mut visit)?;
        }
        Ok(())
    }
}

fn visit_model(
    structure: &Structure,
    model: ModelIndex,
    model_number: i32,
    visit: &mut impl for<'row> FnMut(CanonicalAtomRow<'row>),
) -> Result<(), CifWriteError> {
    for chain in structure.data().chains() {
        let Some(chain_label) = chain.label().filter(|text| !text.is_empty()) else {
            let atom = chain
                .residues()
                .next()
                .and_then(|residue| residue.atoms().next());
            let position = atom.map_or(0, |value| value.index().get());
            return Err(CifWriteError::MissingAtomField {
                atom: position,
                field: "label_asym_id",
            });
        };
        let positions = structure.model_positions(model);
        for residue in chain.residues() {
            for atom in residue.atoms() {
                let row = atom_row(
                    structure,
                    atom,
                    residue,
                    chain,
                    chain_label,
                    positions,
                    model_number,
                )?;
                visit(row);
            }
        }
    }
    Ok(())
}

fn atom_row<'a>(
    structure: &'a Structure,
    atom: AtomRef<'a>,
    residue: ResidueRef<'a>,
    chain: ChainRef<'a>,
    chain_label: &'a str,
    positions: Option<&'a [[f32; 3]]>,
    model_number: i32,
) -> Result<CanonicalAtomRow<'a>, CifWriteError> {
    let atom_index = atom.index().get();
    let Some(atom_site_id) = atom.atom_site_id().filter(|value| *value != 0) else {
        return Err(CifWriteError::MissingAtomSiteId { atom: atom_index });
    };
    let Some(atom_name) = atom.name().filter(|text| !text.is_empty()) else {
        return Err(CifWriteError::MissingAtomField {
            atom: atom_index,
            field: "label_atom_id",
        });
    };
    let Some(component_name) = atom.component_name().filter(|text| !text.is_empty()) else {
        return Err(CifWriteError::MissingAtomField {
            atom: atom_index,
            field: "label_comp_id",
        });
    };
    let position = positions.and_then(|values| values.get(atom.index().as_usize()).copied());
    Ok(CanonicalAtomRow {
        structure,
        atom,
        residue,
        chain,
        atom_site_id,
        atom_name,
        component_name,
        chain_label,
        position,
        model_number,
    })
}

fn model_index(position: usize) -> Result<ModelIndex, CifWriteError> {
    u32::try_from(position)
        .map(ModelIndex::new)
        .map_err(|_| CifWriteError::ModelIndexOverflow { model: position })
}

/// One allocation-free atom-site row over a borrowed model snapshot.
#[derive(Clone, Copy, Debug)]
pub struct CanonicalAtomRow<'a> {
    structure: &'a Structure,
    atom: AtomRef<'a>,
    residue: ResidueRef<'a>,
    chain: ChainRef<'a>,
    atom_site_id: u32,
    atom_name: &'a str,
    component_name: &'a str,
    chain_label: &'a str,
    position: Option<[f32; 3]>,
    model_number: i32,
}

impl<'a> CanonicalAtomRow<'a> {
    /// PDB record family implied by the residue kind.
    #[must_use]
    pub fn group(self) -> &'static str {
        if self.residue.is_het() {
            "HETATM"
        } else {
            "ATOM"
        }
    }

    /// Deposited atom-site identifier.
    #[must_use]
    pub const fn atom_site_id(self) -> i64 {
        self.atom_site_id as i64
    }

    /// Element or an explicit unknown value.
    #[must_use]
    pub fn element(self) -> CanonicalValue<Element> {
        match self.atom.element() {
            Some(value) => CanonicalValue::Present(value),
            None => CanonicalValue::Unknown,
        }
    }

    /// Normalised atom identifier.
    #[must_use]
    pub const fn label_atom_id(self) -> &'a str {
        self.atom_name
    }

    /// Alternate-location identifier, inapplicable for the shared position.
    #[must_use]
    pub fn label_alt_id(self) -> CanonicalValue<&'a str> {
        inapplicable_text(self.atom.alt_label())
    }

    /// Effective component identifier.
    #[must_use]
    pub const fn label_comp_id(self) -> &'a str {
        self.component_name
    }

    /// Normalised chain identifier.
    #[must_use]
    pub const fn label_asym_id(self) -> &'a str {
        self.chain_label
    }

    /// Entity identifier, inapplicable when the chain has none.
    #[must_use]
    pub fn label_entity_id(self) -> CanonicalValue<&'a str> {
        let value = self
            .chain
            .entity()
            .and_then(|entity| self.structure.data().topology.entities.id(entity))
            .and_then(|symbol| self.structure.resolve(symbol));
        inapplicable_text(value)
    }

    /// Canonical sequence position.
    #[must_use]
    pub fn label_seq_id(self) -> CanonicalValue<i64> {
        inapplicable_integer(self.residue.label_seq_id())
    }

    /// Depositor insertion code.
    #[must_use]
    pub fn insertion_code(self) -> CanonicalValue<&'a str> {
        unknown_text(self.residue.ins_code())
    }

    /// Cartesian coordinate rounded to canonical text precision.
    #[must_use]
    pub fn coordinate(self, axis: usize) -> CanonicalValue<f64> {
        match self.position.and_then(|point| point.get(axis).copied()) {
            Some(value) => CanonicalValue::Present(round_decimal(f64::from(value), 1_000.0)),
            None => CanonicalValue::Unknown,
        }
    }

    /// Occupancy rounded to canonical text precision.
    #[must_use]
    pub fn occupancy(self) -> CanonicalValue<f64> {
        optional_measurement(self.atom.occupancy())
    }

    /// Isotropic temperature factor rounded to canonical text precision.
    #[must_use]
    pub fn b_factor(self) -> CanonicalValue<f64> {
        optional_measurement(self.atom.b_factor())
    }

    /// Depositor sequence position.
    #[must_use]
    pub fn auth_seq_id(self) -> CanonicalValue<i64> {
        inapplicable_integer(self.residue.auth_seq_id())
    }

    /// Depositor component identifier.
    #[must_use]
    pub fn auth_comp_id(self) -> CanonicalValue<&'a str> {
        unknown_text(self.residue.auth_name())
    }

    /// Depositor chain identifier.
    #[must_use]
    pub fn auth_asym_id(self) -> CanonicalValue<&'a str> {
        unknown_text(self.chain.auth_label())
    }

    /// Depositor atom identifier.
    #[must_use]
    pub fn auth_atom_id(self) -> CanonicalValue<&'a str> {
        unknown_text(self.atom.auth_name())
    }

    /// Deposited model number.
    #[must_use]
    pub const fn model_number(self) -> i64 {
        self.model_number as i64
    }
}

fn inapplicable_text(value: Option<&str>) -> CanonicalValue<&str> {
    match value.filter(|text| !text.is_empty()) {
        Some(text) => CanonicalValue::Present(text),
        None => CanonicalValue::Inapplicable,
    }
}

fn unknown_text(value: Option<&str>) -> CanonicalValue<&str> {
    match value.filter(|text| !text.is_empty()) {
        Some(text) => CanonicalValue::Present(text),
        None => CanonicalValue::Unknown,
    }
}

fn inapplicable_integer(value: Option<i32>) -> CanonicalValue<i64> {
    match value {
        Some(number) => CanonicalValue::Present(i64::from(number)),
        None => CanonicalValue::Inapplicable,
    }
}

fn optional_measurement(value: Option<f32>) -> CanonicalValue<f64> {
    match value {
        Some(number) => CanonicalValue::Present(round_decimal(f64::from(number), 100.0)),
        None => CanonicalValue::Unknown,
    }
}

pub(crate) fn round_decimal(value: f64, factor: f64) -> f64 {
    (value * factor).round() / factor
}
