//! Stable secondary-structure vocabulary and Python analysis boundary.

use molframe::SecondarySource as Source;
use molframe::SecondaryStructure as Ss;
use pyo3::prelude::*;

/// Declares the Python enum mirroring `Ss` and both conversions from one
/// variant list, so the three can never disagree.
macro_rules! mirror {
    ($($variant:ident),+ $(,)?) => {
        /// Per-residue secondary state; integer values are the native wire codes.
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        #[pyclass(
            name = "SecondaryStructure",
            eq,
            eq_int,
            hash,
            frozen,
            skip_from_py_object,
            module = "molframe"
        )]
        pub(crate) enum PySecondaryStructure {
            $($variant = Ss::$variant.code() as isize,)+
        }

        impl From<Ss> for PySecondaryStructure {
            fn from(value: Ss) -> Self {
                match value {
                    $(Ss::$variant => Self::$variant,)+
                }
            }
        }

        impl PySecondaryStructure {
            const fn native(&self) -> Ss {
                match self {
                    $(Self::$variant => Ss::$variant,)+
                }
            }
        }
    };
}

mirror!(
    Unknown,
    Coil,
    AlphaHelix,
    Strand,
    Turn,
    ThreeTenHelix,
    PiHelix,
    OtherHelix,
    BetaBridge,
    Bend,
    PolyProline,
);

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

/// Where a residue's secondary-structure state came from.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
#[pyclass(
    name = "SecondarySource",
    eq,
    eq_int,
    hash,
    frozen,
    skip_from_py_object,
    module = "molframe"
)]
pub(crate) enum PySecondarySource {
    Unassigned = Source::None.code() as isize,
    File = Source::File.code() as isize,
    Dssp = Source::Dssp.code() as isize,
    CaOnly = Source::CaOnly.code() as isize,
}

impl From<Source> for PySecondarySource {
    fn from(value: Source) -> Self {
        match value {
            Source::File => Self::File,
            Source::Dssp => Self::Dssp,
            Source::CaOnly => Self::CaOnly,
            Source::None => Self::Unassigned,
        }
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
        .map_err(crate::error::failure)?;
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
