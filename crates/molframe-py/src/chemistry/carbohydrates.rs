//! Mechanical conversion of saccharide metadata, rings and linkage reports.

use crate::bindings::{PyStructure, findings_error};
use molframe::chemistry as chem;
use pyo3::prelude::*;
use std::path::PathBuf;

#[derive(Clone, Copy)]
#[pyclass(
    name = "SnfgSymbol",
    module = "molframe.chemistry",
    frozen,
    skip_from_py_object
)]
struct PySnfgSymbol {
    #[pyo3(get)]
    abbreviation: &'static str,
    #[pyo3(get)]
    name: &'static str,
    #[pyo3(get)]
    color: u32,
    #[pyo3(get)]
    shape: &'static str,
    #[pyo3(get)]
    secondary_color: Option<u32>,
}

impl From<chem::SnfgSymbol> for PySnfgSymbol {
    fn from(value: chem::SnfgSymbol) -> Self {
        Self {
            abbreviation: value.abbreviation,
            name: value.name,
            color: value.color,
            shape: value.shape.as_str(),
            secondary_color: value.secondary_color,
        }
    }
}

#[derive(Clone, Copy)]
#[pyclass(
    name = "RingGeometry",
    module = "molframe.chemistry",
    frozen,
    skip_from_py_object
)]
struct PyRingGeometry {
    #[pyo3(get)]
    center: (f32, f32, f32),
    #[pyo3(get)]
    normal: (f32, f32, f32),
    #[pyo3(get)]
    anomeric_direction: Option<(f32, f32, f32)>,
}

impl From<chem::RingGeometry> for PyRingGeometry {
    fn from(value: chem::RingGeometry) -> Self {
        Self {
            center: value.center.into(),
            normal: value.normal.into(),
            anomeric_direction: value.anomeric_direction.map(Into::into),
        }
    }
}

#[derive(Clone)]
#[pyclass(
    name = "Monosaccharide",
    module = "molframe.chemistry",
    frozen,
    skip_from_py_object
)]
struct PyMonosaccharide {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    ring_atoms: Vec<u32>,
    #[pyo3(get)]
    anomeric_atom: Option<u32>,
    #[pyo3(get)]
    symbol: PySnfgSymbol,
    #[pyo3(get)]
    geometry: Option<PyRingGeometry>,
}

impl From<chem::Monosaccharide> for PyMonosaccharide {
    fn from(value: chem::Monosaccharide) -> Self {
        Self {
            residue: value.residue.get(),
            ring_atoms: value
                .ring_atoms
                .into_iter()
                .map(molframe::AtomIndex::get)
                .collect(),
            anomeric_atom: value.anomeric_atom.map(molframe::AtomIndex::get),
            symbol: value.symbol.into(),
            geometry: value.geometry.map(Into::into),
        }
    }
}

#[derive(Clone, Copy)]
#[pyclass(
    name = "CarbohydrateLink",
    module = "molframe.chemistry",
    frozen,
    skip_from_py_object
)]
struct PyCarbohydrateLink {
    #[pyo3(get)]
    donor: usize,
    #[pyo3(get)]
    acceptor: Option<usize>,
    #[pyo3(get)]
    donor_atom: u32,
    #[pyo3(get)]
    acceptor_atom: u32,
    #[pyo3(get)]
    provenance: &'static str,
}

impl From<chem::CarbohydrateLink> for PyCarbohydrateLink {
    fn from(value: chem::CarbohydrateLink) -> Self {
        Self {
            donor: value.donor,
            acceptor: value.acceptor,
            donor_atom: value.donor_atom.get(),
            acceptor_atom: value.acceptor_atom.get(),
            provenance: match value.provenance {
                molframe::BondProvenance::File => "file",
                molframe::BondProvenance::ChemicalComponentDictionary => {
                    "chemical_component_dictionary"
                }
                molframe::BondProvenance::InferredDistance => "inferred_distance",
                molframe::BondProvenance::User => "user",
            },
        }
    }
}

#[pyclass(
    name = "CarbohydrateReport",
    module = "molframe.chemistry",
    frozen,
    skip_from_py_object
)]
struct PyCarbohydrateReport {
    #[pyo3(get)]
    monosaccharides: Vec<PyMonosaccharide>,
    #[pyo3(get)]
    links: Vec<PyCarbohydrateLink>,
    #[pyo3(get)]
    terminal_links: Vec<PyCarbohydrateLink>,
    #[pyo3(get)]
    incomplete_residues: Vec<u32>,
    #[pyo3(get)]
    dictionary_version: Option<String>,
}

impl From<chem::CarbohydrateReport> for PyCarbohydrateReport {
    fn from(value: chem::CarbohydrateReport) -> Self {
        Self {
            monosaccharides: value.monosaccharides.into_iter().map(Into::into).collect(),
            links: value.links.into_iter().map(Into::into).collect(),
            terminal_links: value.terminal_links.into_iter().map(Into::into).collect(),
            incomplete_residues: value
                .incomplete_residues
                .into_iter()
                .map(molframe::ResidueIndex::get)
                .collect(),
            dictionary_version: value
                .dictionary_version
                .map(|version| version.as_str().to_owned()),
        }
    }
}

/// Curated SNFG identity for a common CCD identifier or SNFG abbreviation.
#[pyfunction]
fn snfg_symbol(component: &str) -> Option<PySnfgSymbol> {
    chem::snfg_symbol(component).map(Into::into)
}

/// Perceive sugar rings and glycosidic/protein attachments without mutating input.
#[pyfunction]
#[pyo3(signature = (structure, components=None, *, version="unversioned", spatial_fallback=true))]
fn carbohydrates(
    py: Python<'_>,
    structure: &PyStructure,
    components: Option<PathBuf>,
    version: &str,
    spatial_fallback: bool,
) -> PyResult<PyCarbohydrateReport> {
    let source = structure.inner.clone();
    let version = molframe::DictionaryVersion::new(version);
    let report = py.detach(move || {
        let provider = match components {
            Some(path) => Some(
                molframe::read_component_dictionary(&path, version)
                    .map_err(|findings| findings_error(&findings))?
                    .0,
            ),
            None => None,
        };
        chem::carbohydrates(
            source.engine(),
            provider
                .as_ref()
                .map(|provider| provider as &dyn chem::ComponentProvider),
            chem::CarbohydrateOptions { spatial_fallback },
        )
        .map_err(crate::error::kernel)
    })?;
    Ok(report.into())
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySnfgSymbol>()?;
    module.add_class::<PyRingGeometry>()?;
    module.add_class::<PyMonosaccharide>()?;
    module.add_class::<PyCarbohydrateLink>()?;
    module.add_class::<PyCarbohydrateReport>()?;
    module.add_function(wrap_pyfunction!(snfg_symbol, module)?)?;
    module.add_function(wrap_pyfunction!(carbohydrates, module)?)
}
