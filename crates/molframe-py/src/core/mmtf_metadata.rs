//! Typed MMTF metadata retained on a structure snapshot.

use crate::errors::read_error;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use std::collections::BTreeSet;

#[pyclass(name = "MmtfOptionalField", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyMmtfOptionalField {
    BFactor,
    Occupancy,
    AtomId,
    AltLoc,
    InsCode,
    SequenceIndex,
    ChainName,
    EntityList,
}

impl From<PyMmtfOptionalField> for molframe::pdb::MmtfOptionalField {
    fn from(value: PyMmtfOptionalField) -> Self {
        match value {
            PyMmtfOptionalField::BFactor => Self::BFactor,
            PyMmtfOptionalField::Occupancy => Self::Occupancy,
            PyMmtfOptionalField::AtomId => Self::AtomId,
            PyMmtfOptionalField::AltLoc => Self::AltLoc,
            PyMmtfOptionalField::InsCode => Self::InsCode,
            PyMmtfOptionalField::SequenceIndex => Self::SequenceIndex,
            PyMmtfOptionalField::ChainName => Self::ChainName,
            PyMmtfOptionalField::EntityList => Self::EntityList,
        }
    }
}

impl From<molframe::pdb::MmtfOptionalField> for PyMmtfOptionalField {
    fn from(value: molframe::pdb::MmtfOptionalField) -> Self {
        match value {
            molframe::pdb::MmtfOptionalField::BFactor => Self::BFactor,
            molframe::pdb::MmtfOptionalField::Occupancy => Self::Occupancy,
            molframe::pdb::MmtfOptionalField::AtomId => Self::AtomId,
            molframe::pdb::MmtfOptionalField::AltLoc => Self::AltLoc,
            molframe::pdb::MmtfOptionalField::InsCode => Self::InsCode,
            molframe::pdb::MmtfOptionalField::SequenceIndex => Self::SequenceIndex,
            molframe::pdb::MmtfOptionalField::ChainName => Self::ChainName,
            molframe::pdb::MmtfOptionalField::EntityList => Self::EntityList,
        }
    }
}

#[pyclass(name = "MmtfGroupMetadata", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMmtfGroupMetadata {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    atom_names: Vec<String>,
    #[pyo3(get)]
    elements: Option<Vec<String>>,
    #[pyo3(get)]
    single_letter_code: String,
    #[pyo3(get)]
    chem_comp_type: String,
}

#[pymethods]
impl PyMmtfGroupMetadata {
    #[new]
    fn new(
        name: String,
        atom_names: Vec<String>,
        elements: Option<Vec<String>>,
        single_letter_code: String,
        chem_comp_type: String,
    ) -> Self {
        Self {
            name,
            atom_names,
            elements,
            single_letter_code,
            chem_comp_type,
        }
    }
}

impl From<molframe::pdb::MmtfGroupMetadata> for PyMmtfGroupMetadata {
    fn from(value: molframe::pdb::MmtfGroupMetadata) -> Self {
        Self {
            name: value.name.into(),
            atom_names: value.atom_names.into_iter().map(Into::into).collect(),
            elements: value
                .elements
                .map(|values| values.into_iter().map(Into::into).collect()),
            single_letter_code: value.single_letter_code.into(),
            chem_comp_type: value.chem_comp_type.into(),
        }
    }
}

impl From<PyMmtfGroupMetadata> for molframe::pdb::MmtfGroupMetadata {
    fn from(value: PyMmtfGroupMetadata) -> Self {
        Self {
            name: value.name.into(),
            atom_names: value.atom_names.into_iter().map(Into::into).collect(),
            elements: value
                .elements
                .map(|values| values.into_iter().map(Into::into).collect()),
            single_letter_code: value.single_letter_code.into(),
            chem_comp_type: value.chem_comp_type.into(),
        }
    }
}

#[pyclass(name = "MmtfEntityMetadata", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMmtfEntityMetadata {
    #[pyo3(get)]
    description: String,
    #[pyo3(get)]
    kind: String,
    #[pyo3(get)]
    sequence: String,
}

#[pymethods]
impl PyMmtfEntityMetadata {
    #[new]
    fn new(description: String, kind: String, sequence: String) -> Self {
        Self {
            description,
            kind,
            sequence,
        }
    }
}

impl From<molframe::pdb::MmtfEntityMetadata> for PyMmtfEntityMetadata {
    fn from(value: molframe::pdb::MmtfEntityMetadata) -> Self {
        Self {
            description: value.description.into(),
            kind: value.kind.into(),
            sequence: value.sequence.into(),
        }
    }
}

impl From<PyMmtfEntityMetadata> for molframe::pdb::MmtfEntityMetadata {
    fn from(value: PyMmtfEntityMetadata) -> Self {
        Self {
            description: value.description.into(),
            kind: value.kind.into(),
            sequence: value.sequence.into(),
        }
    }
}

#[pyclass(name = "MmtfMetadata", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMmtfMetadata {
    #[pyo3(get)]
    space_group: Option<String>,
    #[pyo3(get)]
    groups: Vec<PyMmtfGroupMetadata>,
    #[pyo3(get)]
    entities: Vec<PyMmtfEntityMetadata>,
    #[pyo3(get)]
    optional_fields: Vec<PyMmtfOptionalField>,
}

#[pymethods]
impl PyMmtfMetadata {
    #[new]
    #[pyo3(signature = (space_group=None, groups=Vec::new(), entities=Vec::new(), optional_fields=Vec::new()))]
    fn new(
        space_group: Option<String>,
        groups: Vec<PyMmtfGroupMetadata>,
        entities: Vec<PyMmtfEntityMetadata>,
        optional_fields: Vec<PyMmtfOptionalField>,
    ) -> Self {
        Self {
            space_group,
            groups,
            entities,
            optional_fields,
        }
    }
}

impl From<molframe::pdb::MmtfMetadata> for PyMmtfMetadata {
    fn from(value: molframe::pdb::MmtfMetadata) -> Self {
        Self {
            space_group: value.space_group.map(Into::into),
            groups: value.groups.into_iter().map(Into::into).collect(),
            entities: value.entities.into_iter().map(Into::into).collect(),
            optional_fields: value.optional_fields.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<PyMmtfMetadata> for molframe::pdb::MmtfMetadata {
    fn from(value: PyMmtfMetadata) -> Self {
        Self {
            space_group: value.space_group.map(Into::into),
            groups: value.groups.into_iter().map(Into::into).collect(),
            entities: value.entities.into_iter().map(Into::into).collect(),
            optional_fields: value
                .optional_fields
                .into_iter()
                .map(Into::into)
                .collect::<BTreeSet<_>>(),
        }
    }
}

#[pyfunction]
pub(crate) fn mmtf_metadata(py: Python<'_>, structure: &PyStructure) -> Option<PyMmtfMetadata> {
    py.detach(move || -> Option<PyMmtfMetadata> {
        structure
            .structure()
            .extensions()
            .get::<molframe::pdb::MmtfMetadata>(molframe::pdb::MMTF_METADATA_EXTENSION)
            .cloned()
            .map(Into::into)
    })
}

#[pyfunction]
pub(crate) fn with_mmtf_metadata(
    py: Python<'_>,
    structure: &PyStructure,
    metadata: PyMmtfMetadata,
) -> PyStructure {
    py.detach(move || -> PyStructure {
        let metadata: molframe::pdb::MmtfMetadata = metadata.into();
        PyStructure::new(
            structure
                .structure()
                .with_extension(molframe::pdb::MMTF_METADATA_EXTENSION, metadata),
        )
    })
}

#[pyfunction]
pub(crate) fn write_mmtf_with_metadata(
    py: Python<'_>,
    structure: &PyStructure,
    metadata: PyMmtfMetadata,
) -> PyResult<Vec<u8>> {
    let metadata: molframe::pdb::MmtfMetadata = metadata.into();
    let structure = structure
        .structure()
        .with_extension(molframe::pdb::MMTF_METADATA_EXTENSION, metadata);
    py.detach(move || molframe::write_mmtf(&structure))
        .map_err(|findings| read_error(py, &findings))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMmtfOptionalField>()?;
    module.add_class::<PyMmtfGroupMetadata>()?;
    module.add_class::<PyMmtfEntityMetadata>()?;
    module.add_class::<PyMmtfMetadata>()?;
    module.add_function(wrap_pyfunction!(mmtf_metadata, module)?)?;
    module.add_function(wrap_pyfunction!(with_mmtf_metadata, module)?)?;
    module.add_function(wrap_pyfunction!(write_mmtf_with_metadata, module)?)?;
    Ok(())
}
