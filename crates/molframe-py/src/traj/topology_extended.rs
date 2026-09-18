//! Native bindings for reusable HOOMD, LAMMPS and TPR topology records.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;
use std::collections::BTreeMap;

pyo3::create_exception!(_native, HoomdXmlError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, LammpsDataError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, LammpsError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, TprError, pyo3::exceptions::PyValueError);

#[pyclass(name = "HoomdBox", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHoomdBox {
    #[pyo3(get)]
    pub(crate) lengths: [f64; 3],
    #[pyo3(get)]
    pub(crate) tilt: [f64; 3],
}

#[pyclass(name = "HoomdInteraction", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyHoomdInteraction {
    #[pyo3(get)]
    pub(crate) kind: String,
    #[pyo3(get)]
    pub(crate) particles: Vec<u32>,
}

fn hoomd_interaction<const N: usize>(
    value: molframe::traj::HoomdInteraction<N>,
) -> PyHoomdInteraction {
    PyHoomdInteraction {
        kind: value.kind.into(),
        particles: value.particles.into(),
    }
}

#[pyclass(name = "HoomdConfiguration", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyHoomdConfiguration {
    #[pyo3(get)]
    pub(crate) step: Option<u64>,
    #[pyo3(get)]
    pub(crate) dimensions: Option<u8>,
    #[pyo3(get)]
    pub(crate) particle_count: usize,
    #[pyo3(get)]
    pub(crate) cell: Option<PyHoomdBox>,
    #[pyo3(get)]
    pub(crate) positions: Vec<[f64; 3]>,
    #[pyo3(get)]
    pub(crate) images: Vec<[i64; 3]>,
    #[pyo3(get)]
    pub(crate) velocities: Vec<[f64; 3]>,
    #[pyo3(get)]
    pub(crate) accelerations: Vec<[f64; 3]>,
    #[pyo3(get)]
    pub(crate) orientations: Vec<[f64; 4]>,
    #[pyo3(get)]
    pub(crate) types: Vec<String>,
    #[pyo3(get)]
    pub(crate) masses: Vec<f64>,
    #[pyo3(get)]
    pub(crate) charges: Vec<f64>,
    #[pyo3(get)]
    pub(crate) diameters: Vec<f64>,
    #[pyo3(get)]
    pub(crate) bodies: Vec<i64>,
    #[pyo3(get)]
    pub(crate) bonds: Vec<PyHoomdInteraction>,
    #[pyo3(get)]
    pub(crate) angles: Vec<PyHoomdInteraction>,
    #[pyo3(get)]
    pub(crate) dihedrals: Vec<PyHoomdInteraction>,
    #[pyo3(get)]
    pub(crate) impropers: Vec<PyHoomdInteraction>,
    #[pyo3(get)]
    pub(crate) extensions: BTreeMap<String, String>,
}

impl From<molframe::traj::HoomdConfiguration> for PyHoomdConfiguration {
    fn from(value: molframe::traj::HoomdConfiguration) -> Self {
        Self {
            step: value.step,
            dimensions: value.dimensions,
            particle_count: value.particle_count,
            cell: value.cell.map(|cell| PyHoomdBox {
                lengths: cell.lengths,
                tilt: cell.tilt,
            }),
            positions: value.positions,
            images: value.images,
            velocities: value.velocities,
            accelerations: value.accelerations,
            orientations: value.orientations,
            types: value.types.into_iter().map(Into::into).collect(),
            masses: value.masses,
            charges: value.charges,
            diameters: value.diameters,
            bodies: value.bodies,
            bonds: value.bonds.into_iter().map(hoomd_interaction).collect(),
            angles: value.angles.into_iter().map(hoomd_interaction).collect(),
            dihedrals: value.dihedrals.into_iter().map(hoomd_interaction).collect(),
            impropers: value.impropers.into_iter().map(hoomd_interaction).collect(),
            extensions: value
                .extensions
                .into_iter()
                .map(|(key, value)| (key.into(), value.into()))
                .collect(),
        }
    }
}

#[pyfunction]
pub(crate) fn parse_hoomd_xml(py: Python<'_>, source: &str) -> PyResult<PyHoomdConfiguration> {
    py.detach(move || -> PyResult<PyHoomdConfiguration> {
        molframe::traj::parse_hoomd_xml(source)
            .map(Into::into)
            .map_err(|error| HoomdXmlError::new_err(error.to_string()))
    })
}

#[pyclass(name = "LammpsAtomStyle", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyLammpsAtomStyle {
    Atomic,
    Charge,
    Molecular,
    Full,
}

impl TryFrom<molframe::traj::LammpsAtomStyle> for PyLammpsAtomStyle {
    type Error = PyErr;

    fn try_from(value: molframe::traj::LammpsAtomStyle) -> Result<Self, Self::Error> {
        match value {
            molframe::traj::LammpsAtomStyle::Atomic => Ok(Self::Atomic),
            molframe::traj::LammpsAtomStyle::Charge => Ok(Self::Charge),
            molframe::traj::LammpsAtomStyle::Molecular => Ok(Self::Molecular),
            molframe::traj::LammpsAtomStyle::Full => Ok(Self::Full),
            _ => Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "native LAMMPS atom style is newer than this binding version",
            )),
        }
    }
}

#[pyclass(name = "LammpsDataAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLammpsDataAtom {
    #[pyo3(get)]
    pub(crate) id: i64,
    #[pyo3(get)]
    pub(crate) molecule: Option<i64>,
    #[pyo3(get)]
    pub(crate) atom_type: u32,
    #[pyo3(get)]
    pub(crate) charge: Option<f64>,
    #[pyo3(get)]
    pub(crate) position: [f64; 3],
    #[pyo3(get)]
    pub(crate) image: Option<[i32; 3]>,
}

impl From<molframe::traj::LammpsDataAtom> for PyLammpsDataAtom {
    fn from(value: molframe::traj::LammpsDataAtom) -> Self {
        Self {
            id: value.id,
            molecule: value.molecule,
            atom_type: value.atom_type,
            charge: value.charge,
            position: value.position,
            image: value.image,
        }
    }
}

#[pyclass(name = "LammpsInteraction", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLammpsInteraction {
    #[pyo3(get)]
    pub(crate) id: i64,
    #[pyo3(get)]
    pub(crate) interaction_type: u32,
    #[pyo3(get)]
    pub(crate) atoms: Vec<i64>,
}

fn lammps_interaction<const N: usize>(
    value: molframe::traj::LammpsInteraction<N>,
) -> PyLammpsInteraction {
    PyLammpsInteraction {
        id: value.id,
        interaction_type: value.interaction_type,
        atoms: value.atoms.into(),
    }
}

#[pyclass(name = "LammpsDataCell", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyLammpsDataCell {
    #[pyo3(get)]
    pub(crate) lower: [f64; 3],
    #[pyo3(get)]
    pub(crate) upper: [f64; 3],
    #[pyo3(get)]
    pub(crate) tilt: [f64; 3],
}

#[pyclass(name = "LammpsData", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLammpsData {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) atom_style: PyLammpsAtomStyle,
    #[pyo3(get)]
    pub(crate) cell: Option<PyLammpsDataCell>,
    #[pyo3(get)]
    pub(crate) masses: BTreeMap<u32, f64>,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyLammpsDataAtom>,
    #[pyo3(get)]
    pub(crate) bonds: Vec<PyLammpsInteraction>,
    #[pyo3(get)]
    pub(crate) angles: Vec<PyLammpsInteraction>,
    #[pyo3(get)]
    pub(crate) dihedrals: Vec<PyLammpsInteraction>,
    #[pyo3(get)]
    pub(crate) impropers: Vec<PyLammpsInteraction>,
    #[pyo3(get)]
    pub(crate) other_sections: BTreeMap<String, Vec<String>>,
}

impl TryFrom<molframe::traj::LammpsData> for PyLammpsData {
    type Error = PyErr;

    fn try_from(value: molframe::traj::LammpsData) -> Result<Self, Self::Error> {
        Ok(Self {
            title: value.title.into(),
            atom_style: value.atom_style.try_into()?,
            cell: value.cell.map(|cell| PyLammpsDataCell {
                lower: cell.lower,
                upper: cell.upper,
                tilt: cell.tilt,
            }),
            masses: value.masses,
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            bonds: value.bonds.into_iter().map(lammps_interaction).collect(),
            angles: value.angles.into_iter().map(lammps_interaction).collect(),
            dihedrals: value
                .dihedrals
                .into_iter()
                .map(lammps_interaction)
                .collect(),
            impropers: value
                .impropers
                .into_iter()
                .map(lammps_interaction)
                .collect(),
            other_sections: value
                .other_sections
                .into_iter()
                .map(|(key, values)| (key.into(), values.into_iter().map(Into::into).collect()))
                .collect(),
        })
    }
}

#[pyfunction]
pub(crate) fn parse_lammps_data(py: Python<'_>, source: &str) -> PyResult<PyLammpsData> {
    py.detach(move || -> PyResult<PyLammpsData> {
        molframe::traj::parse_lammps_data(source)
            .map_err(|error| LammpsDataError::new_err(error.to_string()))?
            .try_into()
    })
}

#[pyfunction]
pub(crate) fn parse_lammps_dump(py: Python<'_>, source: &str) -> PyResult<Vec<PyTimestep>> {
    py.detach(move || -> PyResult<Vec<PyTimestep>> {
        molframe::traj::parse_lammps_dump(source)
            .map_err(|error| LammpsError::new_err(error.to_string()))?
            .into_iter()
            .map(TryInto::try_into)
            .collect()
    })
}

#[pyclass(name = "TprHeader", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTprHeader {
    #[pyo3(get)]
    pub(crate) producer: String,
    #[pyo3(get)]
    pub(crate) format_version: i32,
    #[pyo3(get)]
    pub(crate) generation: i32,
    #[pyo3(get)]
    pub(crate) precision: u8,
    #[pyo3(get)]
    pub(crate) atom_count: usize,
    #[pyo3(get)]
    pub(crate) tag: String,
    #[pyo3(get)]
    pub(crate) has_coordinates: bool,
}

#[pyclass(name = "TprResidue", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTprResidue {
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) number: i32,
    #[pyo3(get)]
    pub(crate) molecule: usize,
    #[pyo3(get)]
    pub(crate) molecule_type: String,
}

#[pyclass(name = "TprAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTprAtom {
    #[pyo3(get)]
    pub(crate) index: usize,
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) atom_type: String,
    #[pyo3(get)]
    pub(crate) residue: usize,
    #[pyo3(get)]
    pub(crate) mass: f64,
    #[pyo3(get)]
    pub(crate) charge: f64,
    #[pyo3(get)]
    pub(crate) atomic_number: Option<u8>,
}

#[pyclass(name = "TprBond", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyTprBond {
    #[pyo3(get)]
    pub(crate) atom_a: usize,
    #[pyo3(get)]
    pub(crate) atom_b: usize,
}

#[pyclass(name = "TprTopology", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTprTopology {
    #[pyo3(get)]
    pub(crate) header: PyTprHeader,
    #[pyo3(get)]
    pub(crate) residues: Vec<PyTprResidue>,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyTprAtom>,
    #[pyo3(get)]
    pub(crate) bonds: Vec<PyTprBond>,
}

impl From<molframe::traj::TprTopology> for PyTprTopology {
    fn from(value: molframe::traj::TprTopology) -> Self {
        Self {
            header: PyTprHeader {
                producer: value.header.producer,
                format_version: value.header.format_version,
                generation: value.header.generation,
                precision: value.header.precision,
                atom_count: value.header.atom_count,
                tag: value.header.tag,
                has_coordinates: value.header.has_coordinates,
            },
            residues: value
                .residues
                .into_iter()
                .map(|value| PyTprResidue {
                    name: value.name,
                    number: value.number,
                    molecule: value.molecule,
                    molecule_type: value.molecule_type,
                })
                .collect(),
            atoms: value
                .atoms
                .into_iter()
                .map(|value| PyTprAtom {
                    index: value.index,
                    name: value.name,
                    atom_type: value.atom_type,
                    residue: value.residue,
                    mass: value.mass,
                    charge: value.charge,
                    atomic_number: value.atomic_number,
                })
                .collect(),
            bonds: value
                .bonds
                .into_iter()
                .map(|value| PyTprBond {
                    atom_a: value.atom_a,
                    atom_b: value.atom_b,
                })
                .collect(),
        }
    }
}

#[pyfunction]
pub(crate) fn parse_tpr(py: Python<'_>, bytes: Vec<u8>) -> PyResult<PyTprTopology> {
    py.detach(move || -> PyResult<PyTprTopology> {
        molframe::traj::parse_tpr(&bytes)
            .map(Into::into)
            .map_err(|error| TprError::new_err(error.to_string()))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("HoomdXmlError", module.py().get_type::<HoomdXmlError>())?;
    module.add("LammpsDataError", module.py().get_type::<LammpsDataError>())?;
    module.add("LammpsError", module.py().get_type::<LammpsError>())?;
    module.add("TprError", module.py().get_type::<TprError>())?;
    module.add_class::<PyHoomdBox>()?;
    module.add_class::<PyHoomdInteraction>()?;
    module.add_class::<PyHoomdConfiguration>()?;
    module.add_class::<PyLammpsAtomStyle>()?;
    module.add_class::<PyLammpsDataAtom>()?;
    module.add_class::<PyLammpsInteraction>()?;
    module.add_class::<PyLammpsDataCell>()?;
    module.add_class::<PyLammpsData>()?;
    module.add_class::<PyTprHeader>()?;
    module.add_class::<PyTprResidue>()?;
    module.add_class::<PyTprAtom>()?;
    module.add_class::<PyTprBond>()?;
    module.add_class::<PyTprTopology>()?;
    module.add_function(wrap_pyfunction!(parse_hoomd_xml, module)?)?;
    module.add_function(wrap_pyfunction!(parse_lammps_data, module)?)?;
    module.add_function(wrap_pyfunction!(parse_lammps_dump, module)?)?;
    module.add_function(wrap_pyfunction!(parse_tpr, module)?)?;
    Ok(())
}
