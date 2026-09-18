//! Typed access to native chemistry reference data and charge kernels.

#[path = "api.rs"]
mod api;
#[path = "components.rs"]
pub(crate) mod components;
#[path = "conversions.rs"]
mod conversions;
#[path = "equivalence.rs"]
mod equivalence;
#[path = "error.rs"]
mod error;
#[path = "molecule2.rs"]
mod molecule2;
#[path = "molecule_io.rs"]
mod molecule_io;
#[path = "molecules.rs"]
mod molecules;
#[path = "providers.rs"]
mod providers;
#[path = "roles.rs"]
mod roles;
#[path = "side_chain.rs"]
mod side_chain;
#[path = "smarts.rs"]
mod smarts;

pub(crate) use error::value_error;
pub(crate) use roles::PyPolymerRoleProfile;

use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::Arc;

create_exception!(_native, PeoeError, PyValueError);

#[pyclass(name = "ComponentDictionary", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyComponentDictionary(pub(crate) Arc<molframe::CifProvider>);

#[pyfunction]
pub(crate) fn read_component_dictionary(
    py: Python<'_>,
    path: PathBuf,
    version: &str,
) -> PyResult<PyComponentDictionary> {
    // Reading a component dictionary is file-bound work that never touches the
    // interpreter, so other Python threads run while it proceeds.
    py.detach(|| {
        molframe::read_component_dictionary(path, molframe::DictionaryVersion::new(version))
    })
    .map(|(provider, _)| PyComponentDictionary(Arc::new(provider)))
    .map_err(|findings| crate::errors::read_error(py, &findings))
}

#[pymethods]
impl PyComponentDictionary {
    #[new]
    fn new(py: Python<'_>, path: PathBuf, version: &str) -> PyResult<Self> {
        molframe::read_component_dictionary(path, molframe::DictionaryVersion::new(version))
            .map(|(provider, _)| Self(Arc::new(provider)))
            .map_err(|findings| crate::errors::read_error(py, &findings))
    }

    #[getter]
    fn version(&self) -> &str {
        molframe::ComponentProvider::version(self.0.as_ref()).as_str()
    }
}

#[pyclass(name = "Element", frozen, eq, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PyElement(pub(crate) molframe::Element);

#[pymethods]
impl PyElement {
    #[new]
    fn new(symbol: &str) -> PyResult<Self> {
        molframe::Element::from_symbol(symbol)
            .map(Self)
            .ok_or_else(|| PyValueError::new_err(format!("unknown chemical element: {symbol}")))
    }

    #[staticmethod]
    fn from_atomic_number(atomic_number: u8) -> PyResult<Self> {
        let element = molframe::Element::from_atomic_number(atomic_number);
        (!element.is_unknown())
            .then_some(Self(element))
            .ok_or_else(|| PyValueError::new_err("atomic number must be in 1..=118"))
    }

    #[staticmethod]
    fn from_symbol(symbol: &str) -> Option<Self> {
        molframe::Element::from_symbol(symbol).map(Self)
    }

    #[staticmethod]
    fn infer_from_name(name: &str) -> Self {
        Self(molframe::Element::infer_from_name(name))
    }

    #[staticmethod]
    fn infer_from_pdb_atom_name(field: &str) -> Self {
        Self(molframe::Element::infer_from_pdb_atom_name(field))
    }

    #[classattr]
    #[pyo3(name = "UNKNOWN")]
    fn unknown() -> Self {
        Self(molframe::Element::UNKNOWN)
    }

    #[classattr]
    #[pyo3(name = "HYDROGEN")]
    fn hydrogen() -> Self {
        Self(molframe::Element::HYDROGEN)
    }

    #[classattr]
    #[pyo3(name = "CARBON")]
    fn carbon() -> Self {
        Self(molframe::Element::CARBON)
    }

    #[classattr]
    #[pyo3(name = "NITROGEN")]
    fn nitrogen() -> Self {
        Self(molframe::Element::NITROGEN)
    }

    #[classattr]
    #[pyo3(name = "OXYGEN")]
    fn oxygen() -> Self {
        Self(molframe::Element::OXYGEN)
    }

    #[classattr]
    #[pyo3(name = "FLUORINE")]
    fn fluorine() -> Self {
        Self(molframe::Element::FLUORINE)
    }

    #[classattr]
    #[pyo3(name = "PHOSPHORUS")]
    fn phosphorus() -> Self {
        Self(molframe::Element::PHOSPHORUS)
    }

    #[classattr]
    #[pyo3(name = "SULFUR")]
    fn sulfur() -> Self {
        Self(molframe::Element::SULFUR)
    }

    #[classattr]
    #[pyo3(name = "CHLORINE")]
    fn chlorine() -> Self {
        Self(molframe::Element::CHLORINE)
    }

    #[classattr]
    #[pyo3(name = "CALCIUM")]
    fn calcium() -> Self {
        Self(molframe::Element::CALCIUM)
    }

    #[classattr]
    #[pyo3(name = "IRON")]
    fn iron() -> Self {
        Self(molframe::Element::IRON)
    }

    #[classattr]
    #[pyo3(name = "ZINC")]
    fn zinc() -> Self {
        Self(molframe::Element::ZINC)
    }

    #[classattr]
    #[pyo3(name = "SELENIUM")]
    fn selenium() -> Self {
        Self(molframe::Element::SELENIUM)
    }

    #[classattr]
    #[pyo3(name = "MAX_ATOMIC_NUMBER")]
    fn max_atomic_number() -> u8 {
        molframe::Element::MAX_ATOMIC_NUMBER
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
        molframe::element_properties(self.0).map(PyElementProperties::from)
    }

    fn vdw_radius(&self, set: PyRadiusSet) -> Option<f32> {
        molframe::vdw_radius(self.0, set.into())
    }

    fn ionic_radii(&self) -> PyResult<Vec<PyIonicRadius>> {
        molframe::ionic_radii(self.0)
            .map(|values| values.iter().cloned().map(PyIonicRadius::from).collect())
            .map_err(value_error)
    }

    fn is_unknown(&self) -> bool {
        self.0.is_unknown()
    }

    fn is_hydrogen(&self) -> bool {
        self.0.is_hydrogen()
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
    fn table(&self) -> PyRadiusTable {
        let table = molframe::RadiusSet::from(self.clone()).table();
        PyRadiusTable {
            name: table.name.to_owned(),
            version: table.version.to_owned(),
        }
    }

    #[getter]
    fn name(&self) -> &'static str {
        molframe::RadiusSet::from(self.clone()).table().name
    }

    #[getter]
    fn version(&self) -> &'static str {
        molframe::RadiusSet::from(self.clone()).table().version
    }
}

#[pyclass(name = "RadiusTable", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRadiusTable {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    version: String,
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
pub(crate) struct PyPeoeOptions(pub(crate) molframe::PeoeOptions);

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
        Self(molframe::PeoeOptions {
            passes,
            initial_damping,
            damping_factor,
            minimum_electronegativity_difference,
            profile: profile.into(),
        })
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(molframe::PeoeOptions::default())
    }
}

#[pyclass(name = "PeoeParameterProfile", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPeoeParameterProfile {
    GasteigerMarsili,
}

impl From<PyPeoeParameterProfile> for molframe::PeoeParameterProfile {
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

        impl From<PyPeoeAtomType> for molframe::PeoeAtomType {
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
pub(crate) struct PyPeoeAtom(molframe::PeoeAtom);

#[pymethods]
impl PyPeoeAtom {
    #[new]
    fn new(atom_type: PyPeoeAtomType, formal_charge: f64) -> Self {
        Self(molframe::PeoeAtom {
            atom_type: atom_type.into(),
            formal_charge,
        })
    }
}

#[pyclass(name = "PeoeBond", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPeoeBond(molframe::PeoeBond);

#[pymethods]
impl PyPeoeBond {
    #[new]
    fn new(atom_a: usize, atom_b: usize) -> Self {
        Self(molframe::PeoeBond { atom_a, atom_b })
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
    py.detach(move || molframe::peoe_charges(&atoms, &bonds, options.0))
        .map_err(value_error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    components::register(module)?;
    equivalence::register(module)?;
    molecules::register(module)?;
    providers::register(module)?;
    roles::register(module)?;
    side_chain::register(module)?;
    smarts::register(module)?;
    module.add(
        "DEFAULT_BOND_RADIUS_SCALE",
        molframe::DEFAULT_BOND_RADIUS_SCALE,
    )?;
    module.add(
        "DEFAULT_MINIMUM_BOND_DISTANCE",
        molframe::DEFAULT_MINIMUM_BOND_DISTANCE,
    )?;
    module.add_class::<PyComponentDictionary>()?;
    module.add_class::<PyElement>()?;
    module.add_class::<PyRadiusSet>()?;
    module.add("RadiiSet", module.getattr("RadiusSet")?)?;
    module.add_class::<PyRadiusTable>()?;
    module.add_class::<PyElementProperties>()?;
    module.add_class::<PyIonicSpin>()?;
    module.add_class::<PyIonicRadius>()?;
    module.add_class::<PyPeoeOptions>()?;
    module.add_class::<PyPeoeParameterProfile>()?;
    module.add_class::<PyPeoeAtomType>()?;
    module.add_class::<PyPeoeAtom>()?;
    module.add_class::<PyPeoeBond>()?;
    module.add("PeoeError", module.py().get_type::<PeoeError>())?;
    module.add(
        "CifProvider",
        module.py().get_type::<PyComponentDictionary>(),
    )?;
    module.add(
        "ComponentProvider",
        module.py().get_type::<PyComponentDictionary>(),
    )?;
    module.add_function(wrap_pyfunction!(api::component_coverage, module)?)?;
    module.add_function(wrap_pyfunction!(api::component_peoe_charges, module)?)?;
    module.add_function(wrap_pyfunction!(api::apply_component_chemistry, module)?)?;
    module.add_function(wrap_pyfunction!(api::apply_polymer_role_profile, module)?)?;
    module.add_function(wrap_pyfunction!(peoe_charges, module)?)?;
    module.add_function(wrap_pyfunction!(read_component_dictionary, module)?)?;
    Ok(())
}
