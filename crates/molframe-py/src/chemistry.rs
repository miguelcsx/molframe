//! Mechanical adapters for element data and per-atom radii.

use crate::bindings::PyStructure;
use molframe::Element;
use molframe::chemistry::{self as chem, RadiusSet};
use numpy::{PyArray1, ToPyArray};
use pyo3::{exceptions::PyValueError, prelude::*};

pub(crate) fn radius_set(name: &str) -> PyResult<RadiusSet> {
    match name {
        "bondi" => Ok(RadiusSet::Bondi),
        "amber_united" => Ok(RadiusSet::AmberUnited),
        "charmm" => Ok(RadiusSet::Charmm),
        "alvarez" => Ok(RadiusSet::Alvarez),
        _ => Err(PyValueError::new_err(
            "radii must be 'bondi', 'amber_united', 'charmm' or 'alvarez'",
        )),
    }
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
        .ok_or_else(|| PyValueError::new_err(format!("unknown element {symbol:?}")))?;
    let properties = chem::element_properties(element)
        .ok_or_else(|| PyValueError::new_err(format!("no data for element {symbol:?}")))?;
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
        .ok_or_else(|| PyValueError::new_err(format!("unknown element {symbol:?}")))?;
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

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyElementProperties>()?;
    module.add_function(wrap_pyfunction!(element, module)?)?;
    module.add_function(wrap_pyfunction!(vdw_radius, module)?)?;
    module.add_function(wrap_pyfunction!(vdw_radii, module)?)
}
