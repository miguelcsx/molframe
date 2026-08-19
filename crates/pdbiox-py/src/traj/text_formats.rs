//! Native bindings for text trajectory records.

use pyo3::prelude::*;
pyo3::create_exception!(_native, AimsError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, GroError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, TxyzError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, CharmmError, pyo3::exceptions::PyValueError);
#[pyclass(name = "AimsAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAimsAtom {
    #[pyo3(get)]
    pub(crate) species: String,
    #[pyo3(get)]
    pub(crate) position: [f32; 3],
}
#[pymethods]
impl PyAimsAtom {
    #[new]
    fn new(species: String, position: [f32; 3]) -> Self {
        Self { species, position }
    }
}

impl From<PyAimsAtom> for pdbiox::traj::AimsAtom {
    fn from(value: PyAimsAtom) -> Self {
        Self {
            species: value.species.into(),
            position: value.position,
        }
    }
}

impl From<pdbiox::traj::AimsAtom> for PyAimsAtom {
    fn from(value: pdbiox::traj::AimsAtom) -> Self {
        Self {
            species: value.species.into(),
            position: value.position,
        }
    }
}

#[pyclass(name = "AimsGeometry", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAimsGeometry {
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyAimsAtom>,
    #[pyo3(get)]
    pub(crate) lattice_vectors: Option<[[f64; 3]; 3]>,
}

#[pymethods]
impl PyAimsGeometry {
    #[new]
    #[pyo3(signature = (atoms, lattice_vectors=None))]
    fn new(atoms: Vec<PyAimsAtom>, lattice_vectors: Option<[[f64; 3]; 3]>) -> Self {
        Self {
            atoms,
            lattice_vectors,
        }
    }
}

impl From<PyAimsGeometry> for pdbiox::traj::AimsGeometry {
    fn from(value: PyAimsGeometry) -> Self {
        Self {
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            lattice_vectors: value.lattice_vectors,
        }
    }
}

impl From<pdbiox::traj::AimsGeometry> for PyAimsGeometry {
    fn from(value: pdbiox::traj::AimsGeometry) -> Self {
        Self {
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            lattice_vectors: value.lattice_vectors,
        }
    }
}

#[pyclass(name = "GroAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGroAtom {
    #[pyo3(get)]
    pub(crate) residue_number: u32,
    #[pyo3(get)]
    pub(crate) residue_name: String,
    #[pyo3(get)]
    pub(crate) atom_name: String,
    #[pyo3(get)]
    pub(crate) atom_number: u32,
    #[pyo3(get)]
    pub(crate) position: [f32; 3],
    #[pyo3(get)]
    pub(crate) velocity: Option<[f32; 3]>,
}

#[pymethods]
impl PyGroAtom {
    #[new]
    #[pyo3(signature = (residue_number, residue_name, atom_name, atom_number, position, velocity=None))]
    fn new(
        residue_number: u32,
        residue_name: String,
        atom_name: String,
        atom_number: u32,
        position: [f32; 3],
        velocity: Option<[f32; 3]>,
    ) -> Self {
        Self {
            residue_number,
            residue_name,
            atom_name,
            atom_number,
            position,
            velocity,
        }
    }
}

impl From<PyGroAtom> for pdbiox::traj::GroAtom {
    fn from(value: PyGroAtom) -> Self {
        Self {
            residue_number: value.residue_number,
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            atom_number: value.atom_number,
            position: value.position,
            velocity: value.velocity,
        }
    }
}

impl From<pdbiox::traj::GroAtom> for PyGroAtom {
    fn from(value: pdbiox::traj::GroAtom) -> Self {
        Self {
            residue_number: value.residue_number,
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            atom_number: value.atom_number,
            position: value.position,
            velocity: value.velocity,
        }
    }
}

#[pyclass(name = "GroFrame", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGroFrame {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyGroAtom>,
    #[pyo3(get)]
    pub(crate) box_values: Vec<f64>,
}

#[pymethods]
impl PyGroFrame {
    #[new]
    fn new(title: String, atoms: Vec<PyGroAtom>, box_values: Vec<f64>) -> Self {
        Self {
            title,
            atoms,
            box_values,
        }
    }
}

impl From<PyGroFrame> for pdbiox::traj::GroFrame {
    fn from(value: PyGroFrame) -> Self {
        Self {
            title: value.title.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            box_values: value.box_values,
        }
    }
}

impl From<pdbiox::traj::GroFrame> for PyGroFrame {
    fn from(value: pdbiox::traj::GroFrame) -> Self {
        Self {
            title: value.title.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            box_values: value.box_values,
        }
    }
}

#[pyclass(name = "TxyzAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTxyzAtom {
    #[pyo3(get)]
    pub(crate) id: u32,
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) position: [f32; 3],
    #[pyo3(get)]
    pub(crate) atom_type: String,
    #[pyo3(get)]
    pub(crate) bonds: Vec<u32>,
}

#[pymethods]
impl PyTxyzAtom {
    #[new]
    fn new(id: u32, name: String, position: [f32; 3], atom_type: String, bonds: Vec<u32>) -> Self {
        Self {
            id,
            name,
            position,
            atom_type,
            bonds,
        }
    }
}

impl From<PyTxyzAtom> for pdbiox::traj::TxyzAtom {
    fn from(value: PyTxyzAtom) -> Self {
        Self {
            id: value.id,
            name: value.name.into(),
            position: value.position,
            atom_type: value.atom_type.into(),
            bonds: value.bonds,
        }
    }
}

impl From<pdbiox::traj::TxyzAtom> for PyTxyzAtom {
    fn from(value: pdbiox::traj::TxyzAtom) -> Self {
        Self {
            id: value.id,
            name: value.name.into(),
            position: value.position,
            atom_type: value.atom_type.into(),
            bonds: value.bonds,
        }
    }
}

#[pyclass(name = "TxyzFrame", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTxyzFrame {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyTxyzAtom>,
}

#[pymethods]
impl PyTxyzFrame {
    #[new]
    fn new(title: String, atoms: Vec<PyTxyzAtom>) -> Self {
        Self { title, atoms }
    }
}

impl From<PyTxyzFrame> for pdbiox::traj::TxyzFrame {
    fn from(value: PyTxyzFrame) -> Self {
        Self {
            title: value.title.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<pdbiox::traj::TxyzFrame> for PyTxyzFrame {
    fn from(value: pdbiox::traj::TxyzFrame) -> Self {
        Self {
            title: value.title.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyclass(name = "CharmmCardFormat", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCharmmCardFormat {
    Standard,
    Extended,
    Free,
}

impl From<PyCharmmCardFormat> for pdbiox::traj::CharmmCardFormat {
    fn from(value: PyCharmmCardFormat) -> Self {
        match value {
            PyCharmmCardFormat::Standard => Self::Standard,
            PyCharmmCardFormat::Extended => Self::Extended,
            PyCharmmCardFormat::Free => Self::Free,
        }
    }
}

impl From<pdbiox::traj::CharmmCardFormat> for PyCharmmCardFormat {
    fn from(value: pdbiox::traj::CharmmCardFormat) -> Self {
        match value {
            pdbiox::traj::CharmmCardFormat::Standard => Self::Standard,
            pdbiox::traj::CharmmCardFormat::Extended => Self::Extended,
            pdbiox::traj::CharmmCardFormat::Free => Self::Free,
        }
    }
}

#[pyclass(name = "CharmmAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCharmmAtom {
    #[pyo3(get)]
    pub(crate) atom_number: i64,
    #[pyo3(get)]
    pub(crate) residue_number: i64,
    #[pyo3(get)]
    pub(crate) residue_name: String,
    #[pyo3(get)]
    pub(crate) atom_name: String,
    #[pyo3(get)]
    pub(crate) position: [f32; 3],
    #[pyo3(get)]
    pub(crate) segment_id: String,
    #[pyo3(get)]
    pub(crate) residue_id: String,
    #[pyo3(get)]
    pub(crate) weight: f32,
}

#[pymethods]
impl PyCharmmAtom {
    #[new]
    fn new(
        atom_number: i64,
        residue_number: i64,
        residue_name: String,
        atom_name: String,
        position: [f32; 3],
        segment_id: String,
        residue_id: String,
        weight: f32,
    ) -> Self {
        Self {
            atom_number,
            residue_number,
            residue_name,
            atom_name,
            position,
            segment_id,
            residue_id,
            weight,
        }
    }
}

impl From<PyCharmmAtom> for pdbiox::traj::CharmmAtom {
    fn from(value: PyCharmmAtom) -> Self {
        Self {
            atom_number: value.atom_number,
            residue_number: value.residue_number,
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            position: value.position,
            segment_id: value.segment_id.into(),
            residue_id: value.residue_id.into(),
            weight: value.weight,
        }
    }
}

impl From<pdbiox::traj::CharmmAtom> for PyCharmmAtom {
    fn from(value: pdbiox::traj::CharmmAtom) -> Self {
        Self {
            atom_number: value.atom_number,
            residue_number: value.residue_number,
            residue_name: value.residue_name.into(),
            atom_name: value.atom_name.into(),
            position: value.position,
            segment_id: value.segment_id.into(),
            residue_id: value.residue_id.into(),
            weight: value.weight,
        }
    }
}

#[pyclass(name = "CharmmCard", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCharmmCard {
    #[pyo3(get)]
    pub(crate) titles: Vec<String>,
    #[pyo3(get)]
    pub(crate) format: PyCharmmCardFormat,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyCharmmAtom>,
}

#[pymethods]
impl PyCharmmCard {
    #[new]
    fn new(titles: Vec<String>, format: PyCharmmCardFormat, atoms: Vec<PyCharmmAtom>) -> Self {
        Self {
            titles,
            format,
            atoms,
        }
    }
}

impl From<PyCharmmCard> for pdbiox::traj::CharmmCard {
    fn from(value: PyCharmmCard) -> Self {
        Self {
            titles: value.titles.into_iter().map(Into::into).collect(),
            format: value.format.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<pdbiox::traj::CharmmCard> for PyCharmmCard {
    fn from(value: pdbiox::traj::CharmmCard) -> Self {
        Self {
            titles: value.titles.into_iter().map(Into::into).collect(),
            format: value.format.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyfunction]
pub(crate) fn parse_aims_geometry(text: &str) -> PyResult<PyAimsGeometry> {
    pdbiox::traj::parse_aims_geometry(text)
        .map(Into::into)
        .map_err(|error| AimsError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_aims_geometry(geometry: PyAimsGeometry) -> String {
    pdbiox::traj::write_aims_geometry(&geometry.into())
}

#[pyfunction]
pub(crate) fn parse_gro_records(text: &str) -> PyResult<Vec<PyGroFrame>> {
    pdbiox::traj::parse_gro_records(text)
        .map(|frames| frames.into_iter().map(Into::into).collect())
        .map_err(|error| GroError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_gro(frames: Vec<PyGroFrame>) -> PyResult<String> {
    let frames = frames.into_iter().map(Into::into).collect::<Vec<_>>();
    pdbiox::traj::write_gro(&frames).map_err(|error| GroError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn parse_txyz_records(text: &str) -> PyResult<Vec<PyTxyzFrame>> {
    pdbiox::traj::parse_txyz_records(text)
        .map(|frames| frames.into_iter().map(Into::into).collect())
        .map_err(|error| TxyzError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_txyz(frames: Vec<PyTxyzFrame>) -> PyResult<String> {
    let frames = frames.into_iter().map(Into::into).collect::<Vec<_>>();
    pdbiox::traj::write_txyz(&frames).map_err(|error| TxyzError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn parse_charmm_record(text: &str) -> PyResult<PyCharmmCard> {
    pdbiox::traj::parse_charmm_record(text)
        .map(Into::into)
        .map_err(|error| CharmmError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_charmm_card(card: PyCharmmCard) -> PyResult<String> {
    pdbiox::traj::write_charmm_card(&card.into())
        .map_err(|error| CharmmError::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("AimsError", module.py().get_type::<AimsError>())?;
    module.add("GroError", module.py().get_type::<GroError>())?;
    module.add("TxyzError", module.py().get_type::<TxyzError>())?;
    module.add("CharmmError", module.py().get_type::<CharmmError>())?;
    module.add_class::<PyAimsAtom>()?;
    module.add_class::<PyAimsGeometry>()?;
    module.add_class::<PyGroAtom>()?;
    module.add_class::<PyGroFrame>()?;
    module.add_class::<PyTxyzAtom>()?;
    module.add_class::<PyTxyzFrame>()?;
    module.add_class::<PyCharmmCardFormat>()?;
    module.add_class::<PyCharmmAtom>()?;
    module.add_class::<PyCharmmCard>()?;
    module.add_function(wrap_pyfunction!(parse_aims_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(write_aims_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(parse_gro_records, module)?)?;
    module.add_function(wrap_pyfunction!(write_gro, module)?)?;
    module.add_function(wrap_pyfunction!(parse_txyz_records, module)?)?;
    module.add_function(wrap_pyfunction!(write_txyz, module)?)?;
    module.add_function(wrap_pyfunction!(parse_charmm_record, module)?)?;
    module.add_function(wrap_pyfunction!(write_charmm_card, module)?)?;
    Ok(())
}
