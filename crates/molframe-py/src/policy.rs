//! The typed analysis policy: every decision an analysis would otherwise make
//! silently, named once and carried into selections and results.

use crate::bindings::{PySelection, PyStructure};
use molframe::AnalysisPolicy;
use pyo3::prelude::*;

/// An immutable set of analysis decisions.
#[derive(Clone, Debug)]
#[pyclass(name = "AnalysisPolicy", frozen, from_py_object, module = "molframe")]
pub(crate) struct PyAnalysisPolicy(pub(crate) AnalysisPolicy);

/// Parses one decision through the Rust vocabulary that owns its spelling.
fn parse<T>(value: &str) -> PyResult<T>
where
    T: std::str::FromStr<Err = molframe::PolicyParseError>,
{
    value.parse().map_err(crate::error::kernel)
}

/// A decision's canonical spelling in the underscore form Python uses; a
/// `name:payload` value keeps its payload as written.
pub(crate) fn snake(canonical: &str) -> String {
    match canonical.split_once(':') {
        Some((head, payload)) => format!("{}:{payload}", head.replace('-', "_")),
        None => canonical.replace('-', "_"),
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
        assembly=None,
        model=None,
        atom_equivalence=None,
        alignment=None,
        periodic=None,
        contact_def=None,
        float_tolerance=None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        identifiers: Option<&str>,
        altloc: Option<&str>,
        missing_atoms: Option<&str>,
        hydrogens: Option<&str>,
        symmetry: Option<&str>,
        vdw_radii: Option<&str>,
        precision: Option<&str>,
        assembly: Option<&str>,
        model: Option<&str>,
        atom_equivalence: Option<&str>,
        alignment: Option<&str>,
        periodic: Option<&str>,
        contact_def: Option<&str>,
        float_tolerance: Option<(f64, f64)>,
    ) -> PyResult<Self> {
        let mut policy = AnalysisPolicy::default();
        if let Some(value) = assembly {
            policy.assembly = parse(value)?;
        }
        if let Some(value) = model {
            policy.model = parse(value)?;
        }
        if let Some(value) = atom_equivalence {
            policy.atom_equivalence = parse(value)?;
        }
        if let Some(value) = alignment {
            policy.alignment = parse(value)?;
        }
        if let Some(value) = periodic {
            policy.periodic = parse(value)?;
        }
        if let Some(value) = contact_def {
            policy.contact_def = parse(value)?;
        }
        if let Some((relative, absolute)) = float_tolerance {
            if !(relative.is_finite() && absolute.is_finite() && relative >= 0.0 && absolute >= 0.0)
            {
                return Err(crate::error::value(
                    "float_tolerance must be two finite non-negative numbers (relative, absolute)",
                ));
            }
            policy.float_tolerance = molframe::Tolerance { relative, absolute };
        }
        if let Some(value) = identifiers {
            policy.identifiers = parse(value)?;
        }
        if let Some(value) = altloc {
            policy.altloc = parse(value)?;
        }
        if let Some(value) = missing_atoms {
            policy.missing_atoms = parse(value)?;
        }
        if let Some(value) = hydrogens {
            policy.hydrogens = parse(value)?;
        }
        if let Some(value) = symmetry {
            policy.symmetry = parse(value)?;
        }
        if let Some(value) = vdw_radii {
            policy.vdw_radii = parse(value)?;
        }
        if let Some(value) = precision {
            policy.precision = parse(value)?;
        }
        Ok(Self(policy))
    }

    #[getter]
    fn identifiers(&self) -> String {
        snake(self.0.identifiers.name())
    }

    #[getter]
    fn altloc(&self) -> String {
        snake(&self.0.altloc.to_string())
    }

    #[getter]
    fn missing_atoms(&self) -> String {
        snake(self.0.missing_atoms.name())
    }

    #[getter]
    fn hydrogens(&self) -> String {
        snake(self.0.hydrogens.name())
    }

    #[getter]
    fn symmetry(&self) -> String {
        snake(self.0.symmetry.name())
    }

    #[getter]
    fn vdw_radii(&self) -> String {
        snake(self.0.vdw_radii.name())
    }

    #[getter]
    fn precision(&self) -> String {
        snake(self.0.precision.name())
    }

    #[getter]
    fn assembly(&self) -> String {
        snake(&self.0.assembly.to_string())
    }

    #[getter]
    fn model(&self) -> String {
        snake(&self.0.model.to_string())
    }

    #[getter]
    fn atom_equivalence(&self) -> String {
        snake(self.0.atom_equivalence.name())
    }

    #[getter]
    fn alignment(&self) -> String {
        snake(&self.0.alignment.to_string())
    }

    #[getter]
    fn periodic(&self) -> String {
        snake(self.0.periodic.name())
    }

    #[getter]
    fn contact_def(&self) -> String {
        snake(&self.0.contact_def.to_string())
    }

    /// `(relative, absolute)` tolerance for comparing floating-point values.
    #[getter]
    fn float_tolerance(&self) -> (f64, f64) {
        (
            self.0.float_tolerance.relative,
            self.0.float_tolerance.absolute,
        )
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

/// The policy a call names, or the default profile.
pub(crate) fn policy_of(policy: Option<PyRef<'_, PyAnalysisPolicy>>) -> molframe::AnalysisPolicy {
    policy.map_or_else(molframe::AnalysisPolicy::default, |policy| policy.0.clone())
}
