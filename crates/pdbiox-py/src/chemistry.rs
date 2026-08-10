//! Typed access to native chemistry reference data and charge kernels.

use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::Arc;

#[pyclass(name = "ComponentDictionary", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyComponentDictionary(pub(crate) Arc<pdbiox::CifProvider>);

#[pymethods]
impl PyComponentDictionary {
    #[new]
    fn new(py: Python<'_>, path: PathBuf, version: &str) -> PyResult<Self> {
        pdbiox::read_component_dictionary(path, pdbiox::DictionaryVersion::new(version))
            .map(|(provider, _)| Self(Arc::new(provider)))
            .map_err(|findings| crate::errors::read_error(py, &findings))
    }

    #[getter]
    fn version(&self) -> &str {
        pdbiox::ComponentProvider::version(self.0.as_ref()).as_str()
    }

    fn peoe_charges(
        &self,
        py: Python<'_>,
        component_id: &str,
        options: PyPeoeOptions,
    ) -> PyResult<Vec<f64>> {
        let component = pdbiox::ComponentProvider::get(self.0.as_ref(), component_id)
            .map_err(value_error)?
            .ok_or_else(|| PyKeyError::new_err(component_id.to_owned()))?;
        py.detach(move || pdbiox::component_peoe_charges(&component, options.0))
            .map_err(value_error)
    }
}

#[pyclass(name = "Element", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyElement(pdbiox::Element);

#[pymethods]
impl PyElement {
    #[new]
    fn new(symbol: &str) -> PyResult<Self> {
        pdbiox::Element::from_symbol(symbol)
            .map(Self)
            .ok_or_else(|| PyValueError::new_err(format!("unknown chemical element: {symbol}")))
    }

    #[staticmethod]
    fn from_atomic_number(atomic_number: u8) -> PyResult<Self> {
        let element = pdbiox::Element::from_atomic_number(atomic_number);
        (!element.is_unknown())
            .then_some(Self(element))
            .ok_or_else(|| PyValueError::new_err("atomic number must be in 1..=118"))
    }

    #[getter]
    fn symbol(&self) -> &'static str {
        self.0.symbol()
    }

    #[getter]
    const fn atomic_number(&self) -> u8 {
        self.0.atomic_number()
    }

    #[getter]
    fn properties(&self) -> Option<PyElementProperties> {
        pdbiox::element_properties(self.0).map(PyElementProperties::from)
    }

    fn vdw_radius(&self, set: PyRadiusSet) -> Option<f32> {
        pdbiox::vdw_radius(self.0, set.into())
    }

    fn ionic_radii(&self) -> PyResult<Vec<PyIonicRadius>> {
        pdbiox::ionic_radii(self.0)
            .map(|values| values.iter().cloned().map(PyIonicRadius::from).collect())
            .map_err(value_error)
    }

    fn __repr__(&self) -> String {
        format!("Element('{}')", self.0.symbol())
    }
}

#[pyclass(name = "RadiusSet", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PyRadiusSet {
    Bondi,
    AmberUnited,
    Charmm,
    Alvarez,
}

#[pymethods]
impl PyRadiusSet {
    #[getter]
    fn name(&self) -> &'static str {
        pdbiox::RadiusSet::from(self.clone()).table().name
    }

    #[getter]
    fn version(&self) -> &'static str {
        pdbiox::RadiusSet::from(self.clone()).table().version
    }
}

#[pyclass(name = "ElementProperties", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyElementProperties {
    #[pyo3(get)]
    atomic_weight: f64,
    #[pyo3(get)]
    covalent_radius: Option<f32>,
    #[pyo3(get)]
    electronegativity: Option<f32>,
    #[pyo3(get)]
    valence_electrons: u8,
    #[pyo3(get)]
    period: u8,
    #[pyo3(get)]
    group: Option<u8>,
}

#[pyclass(name = "IonicSpin", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyIonicSpin {
    Unspecified,
    High,
    Low,
}

#[pyclass(name = "IonicRadius", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyIonicRadius {
    #[pyo3(get)]
    charge: i8,
    #[pyo3(get)]
    coordination: String,
    #[pyo3(get)]
    spin: PyIonicSpin,
    #[pyo3(get)]
    ionic_radius: f32,
    #[pyo3(get)]
    crystal_radius: f32,
    #[pyo3(get)]
    most_reliable: Option<bool>,
}

#[pyclass(name = "PeoeOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPeoeOptions(pdbiox::PeoeOptions);

#[pymethods]
impl PyPeoeOptions {
    #[new]
    fn new(
        passes: usize,
        initial_damping: f64,
        damping_factor: f64,
        minimum_electronegativity_difference: f64,
        profile: PyPeoeParameterProfile,
    ) -> Self {
        Self(pdbiox::PeoeOptions {
            passes,
            initial_damping,
            damping_factor,
            minimum_electronegativity_difference,
            profile: profile.into(),
        })
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(pdbiox::PeoeOptions::default())
    }
}

#[pyclass(name = "PeoeParameterProfile", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPeoeParameterProfile {
    GasteigerMarsili,
}

impl From<PyPeoeParameterProfile> for pdbiox::PeoeParameterProfile {
    fn from(value: PyPeoeParameterProfile) -> Self {
        match value {
            PyPeoeParameterProfile::GasteigerMarsili => Self::GasteigerMarsili,
        }
    }
}

macro_rules! peoe_atom_types {
    ($($variant:ident),+ $(,)?) => {
        #[pyclass(name = "PeoeAtomType", frozen, eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub(crate) enum PyPeoeAtomType { $($variant),+ }

        impl From<PyPeoeAtomType> for pdbiox::PeoeAtomType {
            fn from(value: PyPeoeAtomType) -> Self {
                match value { $(PyPeoeAtomType::$variant => Self::$variant),+ }
            }
        }
    };
}

peoe_atom_types!(
    H, CSp3, CSp2, CSp, NSp3, NSp2, NSp, OSp3, OSp2, FSp3, ClSp3, BrSp3, ISp3, SSp3, SO, SO2, SSp2,
    PSp3, PSp2, SiSp3, SiSp2, SiSp, BSp3, BSp2, BeSp3, BeSp2, MgSp3, MgSp2, MgSp, AlSp3, AlSp2,
);

#[pyclass(name = "PeoeAtom", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPeoeAtom(pdbiox::PeoeAtom);

#[pymethods]
impl PyPeoeAtom {
    #[new]
    fn new(atom_type: PyPeoeAtomType, formal_charge: f64) -> Self {
        Self(pdbiox::PeoeAtom {
            atom_type: atom_type.into(),
            formal_charge,
        })
    }
}

#[pyclass(name = "PeoeBond", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPeoeBond(pdbiox::PeoeBond);

#[pymethods]
impl PyPeoeBond {
    #[new]
    fn new(atom_a: usize, atom_b: usize) -> Self {
        Self(pdbiox::PeoeBond { atom_a, atom_b })
    }
}

#[pyfunction]
pub(crate) fn peoe_charges(
    py: Python<'_>,
    atoms: Vec<PyPeoeAtom>,
    bonds: Vec<PyPeoeBond>,
    options: PyPeoeOptions,
) -> PyResult<Vec<f64>> {
    let atoms = atoms.into_iter().map(|value| value.0).collect::<Vec<_>>();
    let bonds = bonds.into_iter().map(|value| value.0).collect::<Vec<_>>();
    py.detach(move || pdbiox::peoe_charges(&atoms, &bonds, options.0))
        .map_err(value_error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyComponentDictionary>()?;
    module.add_class::<PyElement>()?;
    module.add_class::<PyRadiusSet>()?;
    module.add_class::<PyElementProperties>()?;
    module.add_class::<PyIonicSpin>()?;
    module.add_class::<PyIonicRadius>()?;
    module.add_class::<PyPeoeOptions>()?;
    module.add_class::<PyPeoeParameterProfile>()?;
    module.add_class::<PyPeoeAtomType>()?;
    module.add_class::<PyPeoeAtom>()?;
    module.add_class::<PyPeoeBond>()?;
    module.add_function(wrap_pyfunction!(peoe_charges, module)?)?;
    Ok(())
}

impl From<PyRadiusSet> for pdbiox::RadiusSet {
    fn from(value: PyRadiusSet) -> Self {
        match value {
            PyRadiusSet::Bondi => Self::Bondi,
            PyRadiusSet::AmberUnited => Self::AmberUnited,
            PyRadiusSet::Charmm => Self::Charmm,
            PyRadiusSet::Alvarez => Self::Alvarez,
        }
    }
}

impl From<pdbiox::ElementProperties> for PyElementProperties {
    fn from(value: pdbiox::ElementProperties) -> Self {
        Self {
            atomic_weight: value.atomic_weight,
            covalent_radius: value.covalent_radius,
            electronegativity: value.electronegativity,
            valence_electrons: value.valence_electrons,
            period: value.period,
            group: value.group,
        }
    }
}

impl From<pdbiox::IonicRadius> for PyIonicRadius {
    fn from(value: pdbiox::IonicRadius) -> Self {
        Self {
            charge: value.charge,
            coordination: value.coordination.into(),
            spin: value.spin.into(),
            ionic_radius: value.ionic_radius,
            crystal_radius: value.crystal_radius,
            most_reliable: value.most_reliable,
        }
    }
}

impl From<pdbiox::IonicSpin> for PyIonicSpin {
    fn from(value: pdbiox::IonicSpin) -> Self {
        match value {
            pdbiox::IonicSpin::Unspecified => Self::Unspecified,
            pdbiox::IonicSpin::High => Self::High,
            pdbiox::IonicSpin::Low => Self::Low,
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
