//! The typed analysis policy: every decision an analysis would otherwise make
//! silently, named once and carried into selections and results.

use crate::bindings::{PySelection, PyStructure};
use molframe::{
    AltlocPolicy, AnalysisPolicy, HydrogenPolicy, MissingPolicy, Namespace, Precision, RadiiSet,
    SymmetryPolicy,
};
use pyo3::{exceptions::PyValueError, prelude::*};

/// An immutable set of analysis decisions.
#[derive(Clone, Debug)]
#[pyclass(name = "AnalysisPolicy", frozen, from_py_object, module = "molframe")]
pub(crate) struct PyAnalysisPolicy(pub(crate) AnalysisPolicy);

fn choose<T: Copy>(name: &str, value: &str, table: &[(&str, T)]) -> PyResult<T> {
    table
        .iter()
        .find_map(|(label, choice)| (*label == value).then_some(*choice))
        .ok_or_else(|| {
            let allowed: Vec<&str> = table.iter().map(|(label, _)| *label).collect();
            PyValueError::new_err(format!("{name} must be one of {}", allowed.join(", ")))
        })
}

fn name_of<T: Copy + PartialEq>(value: T, table: &[(&'static str, T)]) -> &'static str {
    table
        .iter()
        .find_map(|(label, choice)| (*choice == value).then_some(*label))
        .map_or("custom", |label| label)
}

const NAMESPACES: &[(&str, Namespace)] = &[
    ("label", Namespace::Label),
    ("auth", Namespace::Auth),
    ("explicit", Namespace::Explicit),
];
const MISSING: &[(&str, MissingPolicy)] = &[
    ("ignore", MissingPolicy::Ignore),
    ("report", MissingPolicy::Report),
    ("indeterminate", MissingPolicy::Indeterminate),
    ("fail", MissingPolicy::Fail),
];
const HYDROGENS: &[(&str, HydrogenPolicy)] = &[
    ("explicit_only", HydrogenPolicy::ExplicitOnly),
    ("include_inferred", HydrogenPolicy::IncludeInferred),
    ("exclude", HydrogenPolicy::Exclude),
];
const SYMMETRY: &[(&str, SymmetryPolicy)] = &[
    ("none", SymmetryPolicy::None),
    ("crystallographic", SymmetryPolicy::Crystallographic),
    ("biological_assembly", SymmetryPolicy::BiologicalAssembly),
];
const RADII: &[(&str, RadiiSet)] = &[
    ("bondi", RadiiSet::Bondi),
    ("amber_united", RadiiSet::AmberUnited),
    ("charmm", RadiiSet::Charmm),
    ("alvarez", RadiiSet::Alvarez),
];
const PRECISION: &[(&str, Precision)] = &[("f32", Precision::F32), ("f64", Precision::F64)];

fn parse_altloc(value: &str) -> PyResult<AltlocPolicy> {
    if let Some(label) = value.strip_prefix("label:") {
        return if label.is_empty() {
            Err(PyValueError::new_err("altloc 'label:' needs a label"))
        } else {
            Ok(AltlocPolicy::Label(label.into()))
        };
    }
    match value {
        "keep_all" => Ok(AltlocPolicy::KeepAll),
        "conformer_consistent" => Ok(AltlocPolicy::ConformerConsistent),
        "first" => Ok(AltlocPolicy::First),
        "highest_occupancy_per_residue" => Ok(AltlocPolicy::HighestOccupancyPerResidue),
        "highest_occupancy_per_atom" => Ok(AltlocPolicy::HighestOccupancyPerAtom),
        _ => Err(PyValueError::new_err(
            "altloc must be keep_all, conformer_consistent, first, \
             highest_occupancy_per_residue, highest_occupancy_per_atom or label:<id>",
        )),
    }
}

fn altloc_name(policy: &AltlocPolicy) -> String {
    match policy {
        AltlocPolicy::KeepAll => "keep_all".to_owned(),
        AltlocPolicy::ConformerConsistent => "conformer_consistent".to_owned(),
        AltlocPolicy::First => "first".to_owned(),
        AltlocPolicy::HighestOccupancyPerResidue => "highest_occupancy_per_residue".to_owned(),
        AltlocPolicy::HighestOccupancyPerAtom => "highest_occupancy_per_atom".to_owned(),
        AltlocPolicy::Label(label) => format!("label:{label}"),
        _ => "custom".to_owned(),
    }
}

#[pymethods]
impl PyAnalysisPolicy {
    /// Unnamed decisions keep the default profile's choice.
    #[new]
    #[pyo3(signature = (
        *,
        identifiers=None,
        altloc=None,
        missing_atoms=None,
        hydrogens=None,
        symmetry=None,
        vdw_radii=None,
        precision=None,
    ))]
    fn new(
        identifiers: Option<&str>,
        altloc: Option<&str>,
        missing_atoms: Option<&str>,
        hydrogens: Option<&str>,
        symmetry: Option<&str>,
        vdw_radii: Option<&str>,
        precision: Option<&str>,
    ) -> PyResult<Self> {
        let mut policy = AnalysisPolicy::default();
        if let Some(value) = identifiers {
            policy.identifiers = choose("identifiers", value, NAMESPACES)?;
        }
        if let Some(value) = altloc {
            policy.altloc = parse_altloc(value)?;
        }
        if let Some(value) = missing_atoms {
            policy.missing_atoms = choose("missing_atoms", value, MISSING)?;
        }
        if let Some(value) = hydrogens {
            policy.hydrogens = choose("hydrogens", value, HYDROGENS)?;
        }
        if let Some(value) = symmetry {
            policy.symmetry = choose("symmetry", value, SYMMETRY)?;
        }
        if let Some(value) = vdw_radii {
            policy.vdw_radii = choose("vdw_radii", value, RADII)?;
        }
        if let Some(value) = precision {
            policy.precision = choose("precision", value, PRECISION)?;
        }
        Ok(Self(policy))
    }

    #[getter]
    fn identifiers(&self) -> &'static str {
        name_of(self.0.identifiers, NAMESPACES)
    }

    #[getter]
    fn altloc(&self) -> String {
        altloc_name(&self.0.altloc)
    }

    #[getter]
    fn missing_atoms(&self) -> &'static str {
        name_of(self.0.missing_atoms, MISSING)
    }

    #[getter]
    fn hydrogens(&self) -> &'static str {
        name_of(self.0.hydrogens, HYDROGENS)
    }

    #[getter]
    fn symmetry(&self) -> &'static str {
        name_of(self.0.symmetry, SYMMETRY)
    }

    #[getter]
    fn vdw_radii(&self) -> &'static str {
        name_of(self.0.vdw_radii, RADII)
    }

    #[getter]
    fn precision(&self) -> &'static str {
        name_of(self.0.precision, PRECISION)
    }

    /// The named profile when no decision has been changed, else `None`.
    #[getter]
    fn profile(&self) -> Option<&'static str> {
        self.0.profile().map(molframe::ProfileId::as_str)
    }

    /// Stable identity of every decision, for telling results apart.
    #[getter]
    fn fingerprint(&self) -> String {
        self.0.fingerprint().to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    fn __repr__(&self) -> String {
        self.0.profile().map_or_else(
            || format!("AnalysisPolicy(fingerprint={})", self.0.fingerprint()),
            |profile| format!("AnalysisPolicy(profile={profile})"),
        )
    }
}

/// Evaluates a compiled query, raising its warnings and rendering its errors
/// against the query text.
pub(crate) fn select_compiled(
    py: Python<'_>,
    structure: &PyStructure,
    query: &molframe::Query,
    policy: &molframe::AnalysisPolicy,
) -> PyResult<PySelection> {
    use molframe::QueryStructure;

    let evaluation = structure
        .inner
        .select_query(query, policy)
        .map_err(|findings| crate::query_messages::query_error(&findings, query.source()))?;
    crate::query_messages::warn(py, &evaluation.warnings, query.source())?;
    let view = structure.inner.engine().view_of(evaluation.selection);
    Ok(PySelection::from_native(
        structure.clone(),
        molframe::Selection::from(view),
    ))
}
