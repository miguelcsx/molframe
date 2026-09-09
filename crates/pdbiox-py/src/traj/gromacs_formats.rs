//! Native bindings for GROMACS include-topology records.

use pyo3::prelude::*;
use std::collections::BTreeMap;

pyo3::create_exception!(_native, GromacsItpError, pyo3::exceptions::PyValueError);

#[pyclass(name = "GromacsMoleculeType", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGromacsMoleculeType {
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) exclusions: u32,
}

impl From<pdbiox::traj::GromacsMoleculeType> for PyGromacsMoleculeType {
    fn from(value: pdbiox::traj::GromacsMoleculeType) -> Self {
        Self {
            name: value.name.into(),
            exclusions: value.exclusions,
        }
    }
}

impl From<PyGromacsMoleculeType> for pdbiox::traj::GromacsMoleculeType {
    fn from(value: PyGromacsMoleculeType) -> Self {
        Self {
            name: value.name.into(),
            exclusions: value.exclusions,
        }
    }
}

#[pyclass(name = "GromacsItpAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGromacsItpAtom {
    #[pyo3(get)]
    pub(crate) number: u32,
    #[pyo3(get)]
    pub(crate) atom_type: String,
    #[pyo3(get)]
    pub(crate) residue_number: i64,
    #[pyo3(get)]
    pub(crate) residue_name: String,
    #[pyo3(get)]
    pub(crate) atom_name: String,
    #[pyo3(get)]
    pub(crate) charge_group: u32,
    #[pyo3(get)]
    pub(crate) charge: f64,
    #[pyo3(get)]
    pub(crate) mass: Option<f64>,
    #[pyo3(get)]
    pub(crate) state_b: Vec<String>,
}

impl From<pdbiox::traj::GromacsItpAtom> for PyGromacsItpAtom {
    fn from(value: pdbiox::traj::GromacsItpAtom) -> Self {
        Self {
            number: value.number,
            atom_type: value.atom_type.into(),
            residue_number: value.residue_number,
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            charge_group: value.charge_group,
            charge: value.charge,
            mass: value.mass,
            state_b: value.state_b.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<PyGromacsItpAtom> for pdbiox::traj::GromacsItpAtom {
    fn from(value: PyGromacsItpAtom) -> Self {
        Self {
            number: value.number,
            atom_type: value.atom_type.into(),
            residue_number: value.residue_number,
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            charge_group: value.charge_group,
            charge: value.charge,
            mass: value.mass,
            state_b: value.state_b.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyclass(name = "GromacsInteraction", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGromacsInteraction {
    #[pyo3(get)]
    pub(crate) atoms: Vec<u32>,
    #[pyo3(get)]
    pub(crate) function: u32,
    #[pyo3(get)]
    pub(crate) parameters: Vec<String>,
}

impl<const N: usize> From<pdbiox::traj::GromacsInteraction<N>> for PyGromacsInteraction {
    fn from(value: pdbiox::traj::GromacsInteraction<N>) -> Self {
        Self {
            atoms: value.atoms.to_vec(),
            function: value.function,
            parameters: value.parameters.into_iter().map(Into::into).collect(),
        }
    }
}

fn interaction<const N: usize>(
    value: PyGromacsInteraction,
) -> PyResult<pdbiox::traj::GromacsInteraction<N>> {
    let atoms: [u32; N] = value.atoms.try_into().map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(format!(
            "GROMACS interaction requires {N} atom indices"
        ))
    })?;
    Ok(pdbiox::traj::GromacsInteraction {
        atoms,
        function: value.function,
        parameters: value.parameters.into_iter().map(Into::into).collect(),
    })
}

#[pyclass(name = "GromacsItp", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGromacsItp {
    #[pyo3(get)]
    pub(crate) directives: Vec<String>,
    #[pyo3(get)]
    pub(crate) molecule_type: Option<PyGromacsMoleculeType>,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyGromacsItpAtom>,
    #[pyo3(get)]
    pub(crate) bonds: Vec<PyGromacsInteraction>,
    #[pyo3(get)]
    pub(crate) pairs: Vec<PyGromacsInteraction>,
    #[pyo3(get)]
    pub(crate) angles: Vec<PyGromacsInteraction>,
    #[pyo3(get)]
    pub(crate) dihedrals: Vec<PyGromacsInteraction>,
    #[pyo3(get)]
    pub(crate) constraints: Vec<PyGromacsInteraction>,
    #[pyo3(get)]
    pub(crate) exclusions: Vec<Vec<u32>>,
    #[pyo3(get)]
    pub(crate) other_sections: BTreeMap<String, Vec<String>>,
}

impl TryFrom<PyGromacsItp> for pdbiox::traj::GromacsItp {
    type Error = PyErr;

    fn try_from(value: PyGromacsItp) -> Result<Self, Self::Error> {
        Ok(Self {
            directives: value.directives.into_iter().map(Into::into).collect(),
            molecule_type: value.molecule_type.map(Into::into),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            bonds: value
                .bonds
                .into_iter()
                .map(interaction)
                .collect::<PyResult<_>>()?,
            pairs: value
                .pairs
                .into_iter()
                .map(interaction)
                .collect::<PyResult<_>>()?,
            angles: value
                .angles
                .into_iter()
                .map(interaction)
                .collect::<PyResult<_>>()?,
            dihedrals: value
                .dihedrals
                .into_iter()
                .map(interaction)
                .collect::<PyResult<_>>()?,
            constraints: value
                .constraints
                .into_iter()
                .map(interaction)
                .collect::<PyResult<_>>()?,
            exclusions: value.exclusions,
            other_sections: value
                .other_sections
                .into_iter()
                .map(|(key, values)| (key.into(), values.into_iter().map(Into::into).collect()))
                .collect(),
        })
    }
}

impl TryFrom<pdbiox::traj::GromacsItp> for PyGromacsItp {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::GromacsItp) -> Result<Self, Self::Error> {
        Ok(Self {
            directives: value.directives.into_iter().map(Into::into).collect(),
            molecule_type: value.molecule_type.map(Into::into),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            bonds: value.bonds.into_iter().map(Into::into).collect(),
            pairs: value.pairs.into_iter().map(Into::into).collect(),
            angles: value.angles.into_iter().map(Into::into).collect(),
            dihedrals: value.dihedrals.into_iter().map(Into::into).collect(),
            constraints: value.constraints.into_iter().map(Into::into).collect(),
            exclusions: value.exclusions,
            other_sections: value
                .other_sections
                .into_iter()
                .map(|(key, values)| (key.into(), values.into_iter().map(Into::into).collect()))
                .collect(),
        })
    }
}

#[pyfunction]
pub(crate) fn parse_gromacs_itp(py: Python<'_>, text: &str) -> PyResult<PyGromacsItp> {
    py.detach(move || -> PyResult<PyGromacsItp> {
        pdbiox::traj::parse_gromacs_itp(text)
            .map_err(|error| GromacsItpError::new_err(error.to_string()))?
            .try_into()
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("GromacsItpError", module.py().get_type::<GromacsItpError>())?;
    module.add_class::<PyGromacsMoleculeType>()?;
    module.add_class::<PyGromacsItpAtom>()?;
    module.add_class::<PyGromacsInteraction>()?;
    module.add_class::<PyGromacsItp>()?;
    module.add_function(wrap_pyfunction!(parse_gromacs_itp, module)?)?;
    Ok(())
}
