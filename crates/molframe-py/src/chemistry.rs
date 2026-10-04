//! Mechanical adapters for element data and per-atom radii.

mod carbohydrates;
mod charges;
mod roles;
use crate::bindings::{PyStructure, findings_error};
use molframe::Element;
use molframe::chemistry::{self as chem, RadiusSet};
use numpy::{PyArray1, ToPyArray};
use pyo3::prelude::*;
use std::path::PathBuf;

/// The radius set a name selects, through the spelling Rust owns.
pub(crate) fn radius_set(name: &str) -> PyResult<RadiusSet> {
    name.parse().map_err(crate::error::kernel)
}

/// Reference properties of one element.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "ElementProperties",
    frozen,
    skip_from_py_object,
    module = "molframe.chemistry"
)]
struct PyElementProperties {
    symbol: &'static str,
    atomic_number: u8,
    properties: chem::ElementProperties,
}

#[pymethods]
impl PyElementProperties {
    #[getter]
    fn symbol(&self) -> &'static str {
        self.symbol
    }

    #[getter]
    fn atomic_number(&self) -> u8 {
        self.atomic_number
    }

    #[getter]
    fn atomic_weight(&self) -> f64 {
        self.properties.atomic_weight
    }

    #[getter]
    fn covalent_radius(&self) -> Option<f32> {
        self.properties.covalent_radius
    }

    #[getter]
    fn electronegativity(&self) -> Option<f32> {
        self.properties.electronegativity
    }

    #[getter]
    fn valence_electrons(&self) -> u8 {
        self.properties.valence_electrons
    }

    #[getter]
    fn period(&self) -> u8 {
        self.properties.period
    }

    #[getter]
    fn group(&self) -> Option<u8> {
        self.properties.group
    }

    fn __repr__(&self) -> String {
        format!("ElementProperties({})", self.symbol)
    }
}

/// Properties of the element with this symbol, case-insensitively.
#[pyfunction]
fn element(symbol: &str) -> PyResult<PyElementProperties> {
    let element = Element::from_symbol(symbol)
        .ok_or_else(|| crate::error::value(format!("unknown element {symbol:?}")))?;
    let properties = chem::element_properties(element)
        .ok_or_else(|| crate::error::value(format!("no data for element {symbol:?}")))?;
    Ok(PyElementProperties {
        symbol: element.symbol(),
        atomic_number: element.atomic_number(),
        properties,
    })
}

/// Van der Waals radius of one element in ångström, or `None` when the set has none.
#[pyfunction]
#[pyo3(signature = (symbol, *, radii="bondi"))]
fn vdw_radius(symbol: &str, radii: &str) -> PyResult<Option<f32>> {
    let set = radius_set(radii)?;
    let element = Element::from_symbol(symbol)
        .ok_or_else(|| crate::error::value(format!("unknown element {symbol:?}")))?;
    Ok(chem::vdw_radius(element, set))
}

/// Per-atom van der Waals radii; `nan` marks atoms whose element has none.
#[pyfunction]
#[pyo3(signature = (structure, *, radii="bondi"))]
fn vdw_radii<'py>(
    py: Python<'py>,
    structure: &PyStructure,
    radii: &str,
) -> PyResult<Bound<'py, PyArray1<f32>>> {
    let set = radius_set(radii)?;
    let engine = structure.inner.engine();
    let mut values = vec![f32::NAN; engine.atom_count() as usize];
    for atom in engine.data().atoms() {
        let radius = atom
            .element()
            .and_then(|element| chem::vdw_radius(element, set));
        if let (Some(radius), Some(slot)) = (radius, values.get_mut(atom.index().as_usize())) {
            *slot = radius;
        }
    }
    Ok(values.to_pyarray(py))
}

/// Applies a Chemical Component Dictionary to a structure.
///
/// Adds the dictionary's internal bonds and the per-atom roles (charge,
/// hydrogen-bond donor and acceptor, aromaticity, component kind) that the
/// interaction analyses need. File connectivity is kept. Components the
/// dictionary lacks are reported as one warning, not an error.
#[pyfunction]
#[pyo3(signature = (structure, components, *, version="unversioned"))]
fn annotate(
    py: Python<'_>,
    structure: &PyStructure,
    components: PathBuf,
    version: &str,
) -> PyResult<PyStructure> {
    let source = structure.inner.clone();
    let version = molframe::DictionaryVersion::new(version);
    let report = py.detach(move || {
        let (provider, _) = molframe::read_component_dictionary(&components, version)
            .map_err(|findings| findings_error(&findings))?;
        chem::apply_component_chemistry(
            source.engine(),
            &provider,
            chem::PolymerLinkPolicy::Disabled,
        )
        .map_err(crate::error::kernel)
    })?;
    if !report.findings.is_empty() {
        let first: Vec<String> = report
            .findings
            .iter()
            .take(3)
            .map(ToString::to_string)
            .collect();
        PyErr::warn(
            py,
            &py.get_type::<pyo3::exceptions::PyUserWarning>(),
            &std::ffi::CString::new(
                format!(
                    "{} component finding(s), e.g. {}",
                    report.findings.len(),
                    first.join("; ")
                )
                .replace('\0', " "),
            )?,
            1,
        )?;
    }
    Ok(PyStructure::new(molframe::Structure::from(
        report.structure,
    )))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    carbohydrates::register(module)?;
    charges::register(module)?;
    roles::register(module)?;
    module.add_class::<PyElementProperties>()?;
    module.add_function(wrap_pyfunction!(element, module)?)?;
    module.add_function(wrap_pyfunction!(vdw_radius, module)?)?;
    module.add_function(wrap_pyfunction!(vdw_radii, module)?)?;
    module.add_function(wrap_pyfunction!(annotate, module)?)
}
