//! Stable secondary-structure vocabulary and Python analysis boundary.

use molframe::SecondaryStructure as Ss;
use pyo3::prelude::*;

/// Per-residue secondary state; integer values are the native wire codes.
#[derive(Clone, Debug, PartialEq, Eq)]
#[pyclass(
    name = "SecondaryStructure",
    eq,
    eq_int,
    skip_from_py_object,
    module = "molframe"
)]
pub(crate) enum PySecondaryStructure {
    Unknown = Ss::Unknown.code() as isize,
    Coil = Ss::Coil.code() as isize,
    AlphaHelix = Ss::AlphaHelix.code() as isize,
    Strand = Ss::Strand.code() as isize,
    Turn = Ss::Turn.code() as isize,
    ThreeTenHelix = Ss::ThreeTenHelix.code() as isize,
    PiHelix = Ss::PiHelix.code() as isize,
    OtherHelix = Ss::OtherHelix.code() as isize,
    BetaBridge = Ss::BetaBridge.code() as isize,
    Bend = Ss::Bend.code() as isize,
    PolyProline = Ss::PolyProline.code() as isize,
}

impl From<Ss> for PySecondaryStructure {
    fn from(value: Ss) -> Self {
        match value {
            Ss::Unknown => Self::Unknown,
            Ss::Coil => Self::Coil,
            Ss::AlphaHelix => Self::AlphaHelix,
            Ss::Strand => Self::Strand,
            Ss::Turn => Self::Turn,
            Ss::ThreeTenHelix => Self::ThreeTenHelix,
            Ss::PiHelix => Self::PiHelix,
            Ss::OtherHelix => Self::OtherHelix,
            Ss::BetaBridge => Self::BetaBridge,
            Ss::Bend => Self::Bend,
            Ss::PolyProline => Self::PolyProline,
        }
    }
}

impl PySecondaryStructure {
    const fn native(&self) -> Ss {
        match self {
            Self::Unknown => Ss::Unknown,
            Self::Coil => Ss::Coil,
            Self::AlphaHelix => Ss::AlphaHelix,
            Self::Strand => Ss::Strand,
            Self::Turn => Ss::Turn,
            Self::ThreeTenHelix => Ss::ThreeTenHelix,
            Self::PiHelix => Ss::PiHelix,
            Self::OtherHelix => Ss::OtherHelix,
            Self::BetaBridge => Ss::BetaBridge,
            Self::Bend => Ss::Bend,
            Self::PolyProline => Ss::PolyProline,
        }
    }
}

#[pymethods]
impl PySecondaryStructure {
    fn is_helix(&self) -> bool {
        self.native().is_helix()
    }

    fn is_strand(&self) -> bool {
        self.native().is_strand()
    }

    fn is_sheet_like(&self) -> bool {
        self.native().is_sheet_like()
    }
}

#[cfg(feature = "analysis")]
#[pyfunction]
pub(crate) fn dssp(
    py: Python<'_>,
    structure: &crate::bindings::PyStructure,
) -> PyResult<crate::table::PyTable> {
    let structure = structure.inner.clone();
    let rows = py
        .detach(move || {
            molframe::analysis::secondary_structure(
                structure.engine(),
                &molframe::analysis::DsspOptions::default(),
            )
        })
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
    let residues: Vec<_> = rows.residue().iter().map(|index| index.get()).collect();
    let kinds: Vec<_> = rows
        .kind()
        .iter()
        .map(|kind| u32::from(kind.code()))
        .collect();
    Ok(crate::table::TableBuilder::new(py, rows.len())
        .indices("residue", &residues)
        .indices("kind", &kinds)
        .finish())
}
