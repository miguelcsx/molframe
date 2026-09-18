//! Owned Python projections for structure metadata and sequence namespaces.

use crate::analysis::PyMissingResidue;
use crate::hierarchy::PyChain;
use crate::structure::PyStructure;
use molframe::{ChainSequenceExt, EntityKind, PolymerKind};
use pyo3::prelude::*;

#[pyclass(name = "EntityKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyEntityKind {
    Polymer,
    NonPolymer,
    Water,
    Branched,
    Unknown,
}

#[pymethods]
impl PyEntityKind {
    fn __repr__(slf: PyRef<'_, Self>) -> &'static str {
        let value = *slf;
        drop(slf);
        match value {
            Self::Polymer => "EntityKind.Polymer",
            Self::NonPolymer => "EntityKind.NonPolymer",
            Self::Water => "EntityKind.Water",
            Self::Branched => "EntityKind.Branched",
            Self::Unknown => "EntityKind.Unknown",
        }
    }
}

impl PyEntityKind {
    pub(crate) const fn native(&self) -> EntityKind {
        match self {
            Self::Polymer => EntityKind::Polymer,
            Self::NonPolymer => EntityKind::NonPolymer,
            Self::Water => EntityKind::Water,
            Self::Branched => EntityKind::Branched,
            Self::Unknown => EntityKind::Unknown,
        }
    }
}

impl From<EntityKind> for PyEntityKind {
    fn from(value: EntityKind) -> Self {
        match value {
            EntityKind::Polymer => Self::Polymer,
            EntityKind::NonPolymer => Self::NonPolymer,
            EntityKind::Water => Self::Water,
            EntityKind::Branched => Self::Branched,
            _ => Self::Unknown,
        }
    }
}

#[pyclass(name = "PolymerKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPolymerKind {
    None,
    Protein,
    Dna,
    Rna,
    NucleicHybrid,
    Saccharide,
    Other,
}

#[pymethods]
impl PyPolymerKind {
    #[staticmethod]
    fn none() -> Self {
        Self::None
    }

    fn is_polymer(slf: PyRef<'_, Self>) -> bool {
        let value = *slf;
        drop(slf);
        value.native().is_polymer()
    }

    fn is_nucleic(slf: PyRef<'_, Self>) -> bool {
        let value = *slf;
        drop(slf);
        value.native().is_nucleic()
    }
}

impl PyPolymerKind {
    pub(crate) const fn native(self) -> PolymerKind {
        match self {
            Self::None => PolymerKind::None,
            Self::Protein => PolymerKind::Protein,
            Self::Dna => PolymerKind::Dna,
            Self::Rna => PolymerKind::Rna,
            Self::NucleicHybrid => PolymerKind::NucleicHybrid,
            Self::Saccharide => PolymerKind::Saccharide,
            Self::Other => PolymerKind::Other,
        }
    }
}

impl From<PolymerKind> for PyPolymerKind {
    fn from(value: PolymerKind) -> Self {
        match value {
            PolymerKind::None => Self::None,
            PolymerKind::Protein => Self::Protein,
            PolymerKind::Dna => Self::Dna,
            PolymerKind::Rna => Self::Rna,
            PolymerKind::NucleicHybrid => Self::NucleicHybrid,
            PolymerKind::Saccharide => Self::Saccharide,
            _ => Self::Other,
        }
    }
}

#[cfg(test)]
#[path = "metadata_tests.rs"]
mod tests;

#[pyclass(name = "EntryMetadata", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEntryMetadata {
    #[pyo3(get)]
    pub(crate) id: Option<String>,
    #[pyo3(get)]
    pub(crate) title: Option<String>,
    #[pyo3(get)]
    pub(crate) method: Option<String>,
    #[pyo3(get)]
    pub(crate) resolution: Option<f32>,
}

#[pyclass(name = "ReferenceSequence", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReferenceSequence {
    #[pyo3(get)]
    pub(crate) id: String,
    #[pyo3(get)]
    pub(crate) entity_id: String,
    #[pyo3(get)]
    pub(crate) database_name: Option<String>,
    #[pyo3(get)]
    pub(crate) database_code: Option<String>,
    #[pyo3(get)]
    pub(crate) accession: Option<String>,
    #[pyo3(get)]
    pub(crate) one_letter_code: Option<String>,
}

#[pyclass(name = "ReferenceAlignment", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReferenceAlignment {
    #[pyo3(get)]
    pub(crate) id: String,
    #[pyo3(get)]
    pub(crate) reference_id: String,
    #[pyo3(get)]
    pub(crate) chain_ids: Vec<String>,
    #[pyo3(get)]
    pub(crate) canonical: [i32; 2],
    #[pyo3(get)]
    pub(crate) reference: [i32; 2],
}

#[pyclass(name = "SequenceReferences", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySequenceReferences {
    #[pyo3(get)]
    pub(crate) sequences: Vec<PyReferenceSequence>,
    #[pyo3(get)]
    pub(crate) alignments: Vec<PyReferenceAlignment>,
}

#[pyclass(name = "SequenceMapping", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySequenceMapping {
    #[pyo3(get)]
    pub(crate) residue: u32,
    #[pyo3(get)]
    pub(crate) canonical_position: Option<i32>,
    #[pyo3(get)]
    pub(crate) reference_position: Option<i32>,
    #[pyo3(get)]
    pub(crate) reference_id: Option<String>,
}

#[pymethods]
impl PyStructure {
    fn entry_metadata(&self) -> PyEntryMetadata {
        self.structure().data().entry.clone().into()
    }

    fn sequence_references(&self) -> Option<PySequenceReferences> {
        self.structure()
            .data()
            .extensions
            .get::<molframe::SequenceReferences>(molframe::SEQUENCE_REFERENCES_EXTENSION)
            .cloned()
            .map(Into::into)
    }
}

#[pymethods]
impl PyChain {
    fn entity(&self) -> Option<u32> {
        self.inner_chain()
            .and_then(molframe::ChainRef::entity)
            .map(molframe::EntityIndex::get)
    }

    fn entity_kind(&self) -> PyEntityKind {
        let Some(chain) = self.inner_chain() else {
            return PyEntityKind::Unknown;
        };
        chain
            .entity()
            .and_then(|entity| self.inner.data().topology.entities.kind(entity))
            .map_or(PyEntityKind::Unknown, Into::into)
    }

    fn polymer_kind(&self) -> PyPolymerKind {
        self.inner_chain()
            .map_or(PyPolymerKind::None, |chain| chain.polymer_kind().into())
    }

    fn observed_sequence(&self) -> Vec<String> {
        let Some(chain) = self.inner_chain() else {
            return Vec::new();
        };
        chain
            .observed_sequence()
            .into_iter()
            .filter_map(|symbol| self.inner.data().dictionary.resolve(symbol))
            .map(str::to_owned)
            .collect()
    }

    fn canonical_sequence(&self) -> Vec<String> {
        let Some(chain) = self.inner_chain() else {
            return Vec::new();
        };
        chain
            .canonical_sequence()
            .iter()
            .filter_map(|symbol| self.inner.data().dictionary.resolve(*symbol))
            .map(str::to_owned)
            .collect()
    }

    fn reference_sequences(&self) -> Vec<PyReferenceSequence> {
        let Some(chain) = self.inner_chain() else {
            return Vec::new();
        };
        chain
            .reference_sequences()
            .into_iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn sequence_mapping(&self) -> Vec<PySequenceMapping> {
        let Some(chain) = self.inner_chain() else {
            return Vec::new();
        };
        chain
            .sequence_mapping()
            .into_iter()
            .map(Into::into)
            .collect()
    }

    fn missing_residues(&self) -> Vec<PyMissingResidue> {
        let structure = self.inner.data();
        let Some(chain) = self.inner_chain() else {
            return Vec::new();
        };
        chain
            .missing_residues()
            .into_iter()
            .filter_map(|missing| {
                structure
                    .dictionary
                    .resolve(missing.component)
                    .map(|component| PyMissingResidue::from_sequence(missing, component))
            })
            .collect()
    }
}

impl PyChain {
    fn inner_chain(&self) -> Option<molframe::ChainRef<'_>> {
        self.inner.chain(self.index)
    }
}

impl From<molframe::EntryMetadata> for PyEntryMetadata {
    fn from(value: molframe::EntryMetadata) -> Self {
        Self {
            id: value.id.map(str::into_string),
            title: value.title.map(str::into_string),
            method: value.method.map(str::into_string),
            resolution: value.resolution,
        }
    }
}

impl From<molframe::ReferenceSequence> for PyReferenceSequence {
    fn from(value: molframe::ReferenceSequence) -> Self {
        Self {
            id: value.id.into_string(),
            entity_id: value.entity_id.into_string(),
            database_name: value.database_name.map(str::into_string),
            database_code: value.database_code.map(str::into_string),
            accession: value.accession.map(str::into_string),
            one_letter_code: value.one_letter_code.map(str::into_string),
        }
    }
}

impl From<molframe::ReferenceAlignment> for PyReferenceAlignment {
    fn from(value: molframe::ReferenceAlignment) -> Self {
        Self {
            id: value.id.into_string(),
            reference_id: value.reference_id.into_string(),
            chain_ids: value
                .chain_ids
                .iter()
                .map(std::string::ToString::to_string)
                .collect(),
            canonical: value.canonical,
            reference: value.reference,
        }
    }
}

impl From<molframe::SequenceReferences> for PySequenceReferences {
    fn from(value: molframe::SequenceReferences) -> Self {
        Self {
            sequences: value.sequences.into_iter().map(Into::into).collect(),
            alignments: value.alignments.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<molframe::SequenceMapping<'_>> for PySequenceMapping {
    fn from(value: molframe::SequenceMapping<'_>) -> Self {
        Self {
            residue: value.residue.get(),
            canonical_position: value.canonical_position,
            reference_position: value.reference_position,
            reference_id: value.reference_id.map(str::to_owned),
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyEntityKind>()?;
    module.add_class::<PyPolymerKind>()?;
    module.add_class::<PyEntryMetadata>()?;
    module.add_class::<PyReferenceSequence>()?;
    module.add_class::<PyReferenceAlignment>()?;
    module.add_class::<PySequenceReferences>()?;
    module.add_class::<PySequenceMapping>()?;
    module.add(
        "SEQUENCE_REFERENCES_EXTENSION",
        molframe::SEQUENCE_REFERENCES_EXTENSION,
    )?;
    Ok(())
}
