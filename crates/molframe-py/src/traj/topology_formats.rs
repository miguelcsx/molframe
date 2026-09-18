//! Native bindings for AMBER and PSF topology records.

use pyo3::prelude::*;
use std::collections::BTreeMap;

pyo3::create_exception!(_native, AmberTopologyError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, PsfError, pyo3::exceptions::PyValueError);

#[pyclass(name = "AmberSection", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberSection {
    #[pyo3(get)]
    pub(crate) format: String,
    #[pyo3(get)]
    pub(crate) values: Vec<String>,
}

impl From<molframe::traj::AmberSection> for PyAmberSection {
    fn from(value: molframe::traj::AmberSection) -> Self {
        Self {
            format: value.format.into(),
            values: value.values.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<PyAmberSection> for molframe::traj::AmberSection {
    fn from(value: PyAmberSection) -> Self {
        Self {
            format: value.format.into(),
            values: value.values.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyclass(name = "AmberTopologyAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberTopologyAtom {
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) charge: f64,
    #[pyo3(get)]
    pub(crate) mass: f64,
    #[pyo3(get)]
    pub(crate) type_index: u32,
    #[pyo3(get)]
    pub(crate) atom_type: String,
    #[pyo3(get)]
    pub(crate) atomic_number: Option<u8>,
    #[pyo3(get)]
    pub(crate) residue: u32,
}

impl From<molframe::traj::AmberTopologyAtom> for PyAmberTopologyAtom {
    fn from(value: molframe::traj::AmberTopologyAtom) -> Self {
        Self {
            name: value.name.into(),
            charge: value.charge,
            mass: value.mass,
            type_index: value.type_index,
            atom_type: value.atom_type.into(),
            atomic_number: value.atomic_number,
            residue: value.residue,
        }
    }
}

impl From<PyAmberTopologyAtom> for molframe::traj::AmberTopologyAtom {
    fn from(value: PyAmberTopologyAtom) -> Self {
        Self {
            name: value.name.into(),
            charge: value.charge,
            mass: value.mass,
            type_index: value.type_index,
            atom_type: value.atom_type.into(),
            atomic_number: value.atomic_number,
            residue: value.residue,
        }
    }
}

#[pyclass(name = "AmberTopologyResidue", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberTopologyResidue {
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) atom_start: u32,
    #[pyo3(get)]
    pub(crate) atom_end: u32,
}

impl From<molframe::traj::AmberTopologyResidue> for PyAmberTopologyResidue {
    fn from(value: molframe::traj::AmberTopologyResidue) -> Self {
        Self {
            name: value.name.into(),
            atom_start: value.atom_start,
            atom_end: value.atom_end,
        }
    }
}

impl From<PyAmberTopologyResidue> for molframe::traj::AmberTopologyResidue {
    fn from(value: PyAmberTopologyResidue) -> Self {
        Self {
            name: value.name.into(),
            atom_start: value.atom_start,
            atom_end: value.atom_end,
        }
    }
}

macro_rules! topology_tuple {
    ($name:ident, $python_name:literal, $native:ident, $atoms:ty) => {
        #[pyclass(name = $python_name, frozen, from_py_object)]
        #[derive(Clone, Copy, Debug)]
        pub(crate) struct $name {
            #[pyo3(get)]
            pub(crate) atoms: $atoms,
            #[pyo3(get)]
            pub(crate) parameter: u32,
            #[pyo3(get)]
            pub(crate) includes_hydrogen: bool,
        }

        impl From<molframe::traj::$native> for $name {
            fn from(value: molframe::traj::$native) -> Self {
                Self {
                    atoms: value.atoms,
                    parameter: value.parameter,
                    includes_hydrogen: value.includes_hydrogen,
                }
            }
        }

        impl From<$name> for molframe::traj::$native {
            fn from(value: $name) -> Self {
                Self {
                    atoms: value.atoms,
                    parameter: value.parameter,
                    includes_hydrogen: value.includes_hydrogen,
                }
            }
        }
    };
}

topology_tuple!(
    PyAmberTopologyBond,
    "AmberTopologyBond",
    AmberTopologyBond,
    [u32; 2]
);
topology_tuple!(
    PyAmberTopologyAngle,
    "AmberTopologyAngle",
    AmberTopologyAngle,
    [u32; 3]
);

#[pyclass(name = "AmberTopologyDihedral", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAmberTopologyDihedral {
    #[pyo3(get)]
    pub(crate) atoms: [u32; 4],
    #[pyo3(get)]
    pub(crate) parameter: u32,
    #[pyo3(get)]
    pub(crate) improper: bool,
    #[pyo3(get)]
    pub(crate) ignore_end_group: bool,
    #[pyo3(get)]
    pub(crate) includes_hydrogen: bool,
}

impl From<molframe::traj::AmberTopologyDihedral> for PyAmberTopologyDihedral {
    fn from(value: molframe::traj::AmberTopologyDihedral) -> Self {
        Self {
            atoms: value.atoms,
            parameter: value.parameter,
            improper: value.improper,
            ignore_end_group: value.ignore_end_group,
            includes_hydrogen: value.includes_hydrogen,
        }
    }
}

impl From<PyAmberTopologyDihedral> for molframe::traj::AmberTopologyDihedral {
    fn from(value: PyAmberTopologyDihedral) -> Self {
        Self {
            atoms: value.atoms,
            parameter: value.parameter,
            improper: value.improper,
            ignore_end_group: value.ignore_end_group,
            includes_hydrogen: value.includes_hydrogen,
        }
    }
}

#[pyclass(name = "AmberTopology", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberTopology {
    #[pyo3(get)]
    pub(crate) version: Option<String>,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyAmberTopologyAtom>,
    #[pyo3(get)]
    pub(crate) residues: Vec<PyAmberTopologyResidue>,
    #[pyo3(get)]
    pub(crate) bonds: Vec<PyAmberTopologyBond>,
    #[pyo3(get)]
    pub(crate) angles: Vec<PyAmberTopologyAngle>,
    #[pyo3(get)]
    pub(crate) dihedrals: Vec<PyAmberTopologyDihedral>,
    #[pyo3(get)]
    pub(crate) sections: BTreeMap<String, PyAmberSection>,
}

impl From<molframe::traj::AmberTopology> for PyAmberTopology {
    fn from(value: molframe::traj::AmberTopology) -> Self {
        Self {
            version: value.version.map(Into::into),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            residues: value.residues.into_iter().map(Into::into).collect(),
            bonds: value.bonds.into_iter().map(Into::into).collect(),
            angles: value.angles.into_iter().map(Into::into).collect(),
            dihedrals: value.dihedrals.into_iter().map(Into::into).collect(),
            sections: value
                .sections
                .into_iter()
                .map(|(key, value)| (key.into(), value.into()))
                .collect(),
        }
    }
}

#[pyclass(name = "PsfAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPsfAtom {
    #[pyo3(get)]
    pub(crate) segment: String,
    #[pyo3(get)]
    pub(crate) residue_id: String,
    #[pyo3(get)]
    pub(crate) residue_name: String,
    #[pyo3(get)]
    pub(crate) atom_name: String,
    #[pyo3(get)]
    pub(crate) atom_type: String,
    #[pyo3(get)]
    pub(crate) charge: f64,
    #[pyo3(get)]
    pub(crate) mass: f64,
}

impl From<molframe::traj::PsfAtom> for PyPsfAtom {
    fn from(value: molframe::traj::PsfAtom) -> Self {
        Self {
            segment: value.segment.into(),
            residue_id: value.residue_id.into(),
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            atom_type: value.atom_type.into(),
            charge: value.charge,
            mass: value.mass,
        }
    }
}

impl From<PyPsfAtom> for molframe::traj::PsfAtom {
    fn from(value: PyPsfAtom) -> Self {
        Self {
            segment: value.segment.into(),
            residue_id: value.residue_id.into(),
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            atom_type: value.atom_type.into(),
            charge: value.charge,
            mass: value.mass,
        }
    }
}

#[pyclass(name = "PsfTopology", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPsfTopology {
    #[pyo3(get)]
    pub(crate) titles: Vec<String>,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyPsfAtom>,
    #[pyo3(get)]
    pub(crate) bonds: Vec<[u32; 2]>,
    #[pyo3(get)]
    pub(crate) angles: Vec<[u32; 3]>,
    #[pyo3(get)]
    pub(crate) dihedrals: Vec<[u32; 4]>,
    #[pyo3(get)]
    pub(crate) impropers: Vec<[u32; 4]>,
    #[pyo3(get)]
    pub(crate) donors: Vec<[u32; 2]>,
    #[pyo3(get)]
    pub(crate) acceptors: Vec<[u32; 2]>,
    #[pyo3(get)]
    pub(crate) cross_terms: Vec<[u32; 8]>,
}

impl From<molframe::traj::PsfTopology> for PyPsfTopology {
    fn from(value: molframe::traj::PsfTopology) -> Self {
        Self {
            titles: value.titles.into_iter().map(Into::into).collect(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            bonds: value.bonds,
            angles: value.angles,
            dihedrals: value.dihedrals,
            impropers: value.impropers,
            donors: value.donors,
            acceptors: value.acceptors,
            cross_terms: value.cross_terms,
        }
    }
}

impl From<PyPsfTopology> for molframe::traj::PsfTopology {
    fn from(value: PyPsfTopology) -> Self {
        Self {
            titles: value.titles.into_iter().map(Into::into).collect(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            bonds: value.bonds,
            angles: value.angles,
            dihedrals: value.dihedrals,
            impropers: value.impropers,
            donors: value.donors,
            acceptors: value.acceptors,
            cross_terms: value.cross_terms,
        }
    }
}

#[pyfunction]
pub(crate) fn parse_amber_topology(py: Python<'_>, text: &str) -> PyResult<PyAmberTopology> {
    py.detach(move || -> PyResult<PyAmberTopology> {
        molframe::traj::parse_amber_topology(text)
            .map(Into::into)
            .map_err(|error| AmberTopologyError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn parse_psf(py: Python<'_>, text: &str) -> PyResult<PyPsfTopology> {
    py.detach(move || -> PyResult<PyPsfTopology> {
        molframe::traj::parse_psf(text)
            .map(Into::into)
            .map_err(|error| PsfError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_psf(py: Python<'_>, topology: PyPsfTopology) -> String {
    py.detach(move || -> String { molframe::traj::write_psf(&topology.into()) })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "AmberTopologyError",
        module.py().get_type::<AmberTopologyError>(),
    )?;
    module.add("PsfError", module.py().get_type::<PsfError>())?;
    module.add_class::<PyAmberSection>()?;
    module.add_class::<PyAmberTopologyAtom>()?;
    module.add_class::<PyAmberTopologyResidue>()?;
    module.add_class::<PyAmberTopologyBond>()?;
    module.add_class::<PyAmberTopologyAngle>()?;
    module.add_class::<PyAmberTopologyDihedral>()?;
    module.add_class::<PyAmberTopology>()?;
    module.add_class::<PyPsfAtom>()?;
    module.add_class::<PyPsfTopology>()?;
    module.add_function(wrap_pyfunction!(parse_amber_topology, module)?)?;
    module.add_function(wrap_pyfunction!(parse_psf, module)?)?;
    module.add_function(wrap_pyfunction!(write_psf, module)?)?;
    Ok(())
}
