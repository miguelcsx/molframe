//! Declarative motif and constraint values.

use super::errors::motif_error;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "AtomSite", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAtomSite(pub(crate) pdbiox::fx::AtomSite);

#[pymethods]
impl PyAtomSite {
    #[new]
    fn new(component: String, atom: String) -> PyResult<Self> {
        if component.is_empty() || atom.is_empty() {
            return Err(PyValueError::new_err(
                "component and atom names must not be empty",
            ));
        }
        Ok(Self(pdbiox::fx::AtomSite::new(component, atom)))
    }

    #[getter]
    fn component(&self) -> &str {
        &self.0.component
    }

    #[getter]
    fn atom(&self) -> &str {
        &self.0.atom
    }

    fn __repr__(&self) -> String {
        format!("AtomSite({:?}, {:?})", self.0.component, self.0.atom)
    }

    fn __str__(&self) -> String {
        format!("{}.{}", self.0.component, self.0.atom)
    }
}

#[pyclass(name = "ComponentRole", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyComponentRole {
    Residue,
    Ligand,
    Cofactor,
    Metal,
    Unknown,
}

impl From<PyComponentRole> for pdbiox::fx::ComponentRole {
    fn from(value: PyComponentRole) -> Self {
        match value {
            PyComponentRole::Residue | PyComponentRole::Unknown => Self::Residue,
            PyComponentRole::Ligand => Self::Ligand,
            PyComponentRole::Cofactor => Self::Cofactor,
            PyComponentRole::Metal => Self::Metal,
        }
    }
}

impl From<pdbiox::fx::ComponentRole> for PyComponentRole {
    fn from(value: pdbiox::fx::ComponentRole) -> Self {
        match value {
            pdbiox::fx::ComponentRole::Residue => Self::Residue,
            pdbiox::fx::ComponentRole::Ligand => Self::Ligand,
            pdbiox::fx::ComponentRole::Cofactor => Self::Cofactor,
            pdbiox::fx::ComponentRole::Metal => Self::Metal,
            _ => Self::Unknown,
        }
    }
}

#[pyclass(name = "ComponentSpec", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyComponentSpec(pub(crate) pdbiox::fx::ComponentSpec);

#[pymethods]
impl PyComponentSpec {
    #[new]
    fn new(role: PyComponentRole) -> Self {
        Self(pdbiox::fx::ComponentSpec::new(role.into()))
    }

    #[staticmethod]
    fn residue() -> Self {
        Self::new(PyComponentRole::Residue)
    }

    #[staticmethod]
    fn ligand() -> Self {
        Self::new(PyComponentRole::Ligand)
    }

    #[staticmethod]
    fn cofactor() -> Self {
        Self::new(PyComponentRole::Cofactor)
    }

    #[staticmethod]
    fn metal() -> Self {
        Self::new(PyComponentRole::Metal)
    }

    #[pyo3(signature = (id))]
    fn component(&self, id: String) -> Self {
        Self(self.0.clone().component(id))
    }

    fn require(&self, atom: String) -> Self {
        Self(self.0.clone().require(atom))
    }

    fn equivalent(&self, atoms: Vec<String>) -> Self {
        Self(self.0.clone().equivalent(atoms))
    }

    #[getter]
    fn role(&self) -> PyComponentRole {
        self.0.role.into()
    }

    #[getter]
    fn component_ids(&self) -> Vec<String> {
        self.0
            .component_ids
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    #[getter]
    fn required_atoms(&self) -> Vec<String> {
        self.0
            .required_atoms
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    #[getter]
    fn equivalent_atoms(&self) -> Vec<Vec<String>> {
        self.0
            .equivalent_atoms
            .iter()
            .map(|group| group.iter().map(ToString::to_string).collect())
            .collect()
    }
}

#[pyclass(name = "Constraint", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyConstraint(pub(crate) pdbiox::fx::Constraint);

#[pymethods]
impl PyConstraint {
    #[staticmethod]
    fn distance(
        first: &PyAtomSite,
        second: &PyAtomSite,
        target: f64,
        tolerance: f64,
    ) -> PyResult<Self> {
        validate_tolerance(target, tolerance)?;
        Ok(Self(pdbiox::fx::Constraint::Distance {
            first: first.0.clone(),
            second: second.0.clone(),
            target,
            tolerance,
        }))
    }

    #[staticmethod]
    fn angle(atoms: Vec<PyAtomSite>, target: f64, tolerance: f64) -> PyResult<Self> {
        validate_tolerance(target, tolerance)?;
        let atoms = fixed_sites::<3>(atoms, "angle")?;
        Ok(Self(pdbiox::fx::Constraint::Angle {
            atoms,
            target,
            tolerance,
        }))
    }

    #[staticmethod]
    fn dihedral(atoms: Vec<PyAtomSite>, target: f64, tolerance: f64) -> PyResult<Self> {
        validate_tolerance(target, tolerance)?;
        let atoms = fixed_sites::<4>(atoms, "dihedral")?;
        Ok(Self(pdbiox::fx::Constraint::Dihedral {
            atoms,
            target,
            tolerance,
        }))
    }

    #[staticmethod]
    fn chirality(atoms: Vec<PyAtomSite>, positive: bool) -> PyResult<Self> {
        let atoms = fixed_sites::<4>(atoms, "chirality")?;
        Ok(Self(pdbiox::fx::Constraint::Chirality { atoms, positive }))
    }

    #[staticmethod]
    fn planarity(atoms: Vec<PyAtomSite>, tolerance: f64) -> PyResult<Self> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(PyValueError::new_err(
                "planarity tolerance must be finite and non-negative",
            ));
        }
        if atoms.len() < 3 {
            return Err(PyValueError::new_err(
                "planarity requires at least three atom sites",
            ));
        }
        Ok(Self(pdbiox::fx::Constraint::Planarity {
            atoms: atoms.into_iter().map(|site| site.0).collect(),
            tolerance,
        }))
    }

    #[staticmethod]
    fn coordination(
        centre: &PyAtomSite,
        partners: Vec<PyAtomSite>,
        count: usize,
        max_distance: f64,
    ) -> PyResult<Self> {
        if !max_distance.is_finite() || max_distance < 0.0 {
            return Err(PyValueError::new_err(
                "max_distance must be finite and non-negative",
            ));
        }
        if count > partners.len() {
            return Err(PyValueError::new_err(
                "coordination count cannot exceed the number of partners",
            ));
        }
        Ok(Self(pdbiox::fx::Constraint::Coordination {
            centre: centre.0.clone(),
            partners: partners.into_iter().map(|site| site.0).collect(),
            count,
            max_distance,
        }))
    }

    #[staticmethod]
    fn steric_exclusion(
        first: &PyAtomSite,
        second: &PyAtomSite,
        min_distance: f64,
    ) -> PyResult<Self> {
        if !min_distance.is_finite() || min_distance < 0.0 {
            return Err(PyValueError::new_err(
                "min_distance must be finite and non-negative",
            ));
        }
        Ok(Self(pdbiox::fx::Constraint::StericExclusion {
            first: first.0.clone(),
            second: second.0.clone(),
            min_distance,
        }))
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match &self.0 {
            pdbiox::fx::Constraint::Distance { .. } => "distance",
            pdbiox::fx::Constraint::Angle { .. } => "angle",
            pdbiox::fx::Constraint::Dihedral { .. } => "dihedral",
            pdbiox::fx::Constraint::Chirality { .. } => "chirality",
            pdbiox::fx::Constraint::Planarity { .. } => "planarity",
            pdbiox::fx::Constraint::Coordination { .. } => "coordination",
            pdbiox::fx::Constraint::StericExclusion { .. } => "steric-exclusion",
            _ => "unknown",
        }
    }

    #[getter]
    fn atoms(&self) -> Vec<PyAtomSite> {
        if let pdbiox::fx::Constraint::Distance { first, second, .. }
        | pdbiox::fx::Constraint::StericExclusion { first, second, .. } = &self.0
        {
            return vec![PyAtomSite(first.clone()), PyAtomSite(second.clone())];
        }
        if let pdbiox::fx::Constraint::Angle { atoms, .. } = &self.0 {
            return atoms.iter().cloned().map(PyAtomSite).collect();
        }
        if let pdbiox::fx::Constraint::Dihedral { atoms, .. }
        | pdbiox::fx::Constraint::Chirality { atoms, .. } = &self.0
        {
            return atoms.iter().cloned().map(PyAtomSite).collect();
        }
        if let pdbiox::fx::Constraint::Planarity { atoms, .. } = &self.0 {
            return atoms.iter().cloned().map(PyAtomSite).collect();
        }
        if let pdbiox::fx::Constraint::Coordination {
            centre, partners, ..
        } = &self.0
        {
            return std::iter::once(centre)
                .chain(partners)
                .cloned()
                .map(PyAtomSite)
                .collect();
        }
        Vec::new()
    }

    #[getter]
    fn target(&self) -> Option<f64> {
        match &self.0 {
            pdbiox::fx::Constraint::Distance { target, .. }
            | pdbiox::fx::Constraint::Angle { target, .. }
            | pdbiox::fx::Constraint::Dihedral { target, .. } => Some(*target),
            _ => None,
        }
    }

    #[getter]
    fn tolerance(&self) -> Option<f64> {
        match &self.0 {
            pdbiox::fx::Constraint::Distance { tolerance, .. }
            | pdbiox::fx::Constraint::Angle { tolerance, .. }
            | pdbiox::fx::Constraint::Dihedral { tolerance, .. }
            | pdbiox::fx::Constraint::Planarity { tolerance, .. } => Some(*tolerance),
            _ => None,
        }
    }

    #[getter]
    fn positive(&self) -> Option<bool> {
        match &self.0 {
            pdbiox::fx::Constraint::Chirality { positive, .. } => Some(*positive),
            _ => None,
        }
    }

    #[getter]
    fn count(&self) -> Option<usize> {
        match &self.0 {
            pdbiox::fx::Constraint::Coordination { count, .. } => Some(*count),
            _ => None,
        }
    }

    #[getter]
    fn max_distance(&self) -> Option<f64> {
        match &self.0 {
            pdbiox::fx::Constraint::Coordination { max_distance, .. } => Some(*max_distance),
            _ => None,
        }
    }

    #[getter]
    fn min_distance(&self) -> Option<f64> {
        match &self.0 {
            pdbiox::fx::Constraint::StericExclusion { min_distance, .. } => Some(*min_distance),
            _ => None,
        }
    }
}

#[pyclass(name = "NamedConstraint", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNamedConstraint(pub(crate) pdbiox::fx::NamedConstraint);

#[pymethods]
impl PyNamedConstraint {
    #[new]
    fn new(name: String, constraint: &PyConstraint) -> Self {
        Self(pdbiox::fx::NamedConstraint {
            name: name.into_boxed_str(),
            constraint: constraint.0.clone(),
        })
    }

    #[getter]
    fn name(&self) -> &str {
        &self.0.name
    }

    #[getter]
    fn constraint(&self) -> PyConstraint {
        PyConstraint(self.0.constraint.clone())
    }
}

#[pyclass(name = "Motif", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMotif(pub(crate) pdbiox::fx::Motif);

#[pymethods]
impl PyMotif {
    #[new]
    fn new(
        components: BTreeMap<String, PyComponentSpec>,
        constraints: Vec<PyNamedConstraint>,
    ) -> PyResult<Self> {
        let components = components
            .into_iter()
            .map(|(name, spec)| (name.into_boxed_str(), spec.0))
            .collect::<Vec<_>>();
        let constraints = constraints
            .into_iter()
            .map(|item| item.0)
            .collect::<Vec<_>>();
        pdbiox::fx::Motif::new(components, constraints)
            .map(Self)
            .map_err(motif_error)
    }

    #[getter]
    fn components(&self) -> BTreeMap<String, PyComponentSpec> {
        self.0
            .components()
            .iter()
            .map(|(name, spec)| (name.to_string(), PyComponentSpec(spec.clone())))
            .collect()
    }

    #[getter]
    fn constraints(&self) -> Vec<PyNamedConstraint> {
        self.0
            .constraints()
            .iter()
            .cloned()
            .map(PyNamedConstraint)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "Motif(components={}, constraints={})",
            self.0.components().len(),
            self.0.constraints().len()
        )
    }
}

fn validate_tolerance(target: f64, tolerance: f64) -> PyResult<()> {
    if !target.is_finite() || !tolerance.is_finite() || tolerance < 0.0 {
        return Err(PyValueError::new_err(
            "target must be finite and tolerance must be finite and non-negative",
        ));
    }
    Ok(())
}

fn fixed_sites<const N: usize>(
    values: Vec<PyAtomSite>,
    kind: &str,
) -> PyResult<[pdbiox::fx::AtomSite; N]> {
    if values.len() != N {
        return Err(PyValueError::new_err(format!(
            "{kind} requires exactly {N} atom sites"
        )));
    }
    values
        .into_iter()
        .map(|site| site.0)
        .collect::<Vec<_>>()
        .try_into()
        .map_err(|_| PyValueError::new_err(format!("{kind} received the wrong number of sites")))
}
