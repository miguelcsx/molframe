//! Mechanical adapters for element data and per-atom radii.

mod carbohydrates;
mod charges;
mod linkage;
mod molecules;
mod roles;
use crate::bindings::{PyStructure, findings_error};
use molframe::Element;
use molframe::chemistry::{self as chem, RadiusSet};
use numpy::{PyArray1, ToPyArray};
use pyo3::prelude::*;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::sync::Arc;

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
    let source = structure.inner.clone();
    let values = py.detach(move || chem::atom_radii(source.engine(), set));
    Ok(values.to_pyarray(py))
}

/// Applies a Chemical Component Dictionary to a structure.
///
/// Adds the dictionary's internal bonds and the per-atom roles (charge,
/// hydrogen-bond donor and acceptor, aromaticity, component kind) that the
/// interaction analyses need. File connectivity is kept. Components the
/// dictionary lacks are reported as one warning, not an error.
/// Inter-residue bonds are added only with an explicit `polymer_links` policy.
///
/// Components resolve from the structure's own `chem_comp_*` rows first, then
/// the `components` dictionary when given, then the compact bundled core set
/// (standard and common modified residues, nucleotides, water, common ions,
/// frequent cofactors, carbohydrates); every other component is reported as
/// unknown.
#[pyfunction]
#[pyo3(signature = (
    structure,
    components=None,
    *,
    version="unversioned",
    polymer_links=None,
))]
fn annotate(
    py: Python<'_>,
    structure: &PyStructure,
    components: Option<PathBuf>,
    version: &str,
    polymer_links: Option<PyRef<'_, linkage::PyPolymerLinkPolicy>>,
) -> PyResult<PyStructure> {
    let source = structure.inner.clone();
    let policy = match polymer_links {
        Some(policy) => policy.inner.clone(),
        None => chem::PolymerLinkPolicy::Disabled,
    };
    let version = molframe::DictionaryVersion::new(version);
    let report = py.detach(move || {
        // Only the components the structure uses are parsed and lowered.
        let user: Option<Arc<dyn chem::ComponentProvider>> = match components {
            None => None,
            Some(path) => Some(Arc::new(
                molframe::open_component_dictionary(&path, version)
                    .map_err(|findings| findings_error(&findings))?,
            )),
        };
        let provider = chem::resolution_stack(source.engine(), user);
        chem::apply_component_chemistry(source.engine(), &provider, policy)
            .map_err(crate::error::kernel)
    })?;
    if !report.findings.is_empty() {
        let unknown: Vec<&str> = report
            .findings
            .iter()
            .filter(|finding| finding.code() == molframe::Code::W3201)
            .flat_map(molframe::Diagnostic::context)
            .filter(|item| item.label() == "component")
            .map(molframe::ContextItem::value)
            .collect();
        let first: Vec<String> = report
            .findings
            .iter()
            .filter(|finding| finding.code() != molframe::Code::W3201)
            .take(3)
            .map(ToString::to_string)
            .collect();
        let mut message = format!("{} component finding(s)", report.findings.len());
        if !unknown.is_empty() {
            let _ = write!(
                message,
                "; no definition for component(s) {}: pass the full wwPDB dictionary \
                 (components.cif.gz) as `components` or embed chem_comp_atom and \
                 chem_comp_bond rows in the file",
                unknown.join(", ")
            );
        }
        if !first.is_empty() {
            let _ = write!(message, "; e.g. {}", first.join("; "));
        }
        PyErr::warn(
            py,
            &py.get_type::<pyo3::exceptions::PyUserWarning>(),
            &std::ffi::CString::new(message.replace('\0', " "))?,
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
    linkage::register(module)?;
    molecules::register(module)?;
    roles::register(module)?;
    module.add_class::<PyElementProperties>()?;
    module.add_function(wrap_pyfunction!(element, module)?)?;
    module.add_function(wrap_pyfunction!(vdw_radius, module)?)?;
    module.add_function(wrap_pyfunction!(vdw_radii, module)?)?;
    module.add_function(wrap_pyfunction!(annotate, module)?)
}
