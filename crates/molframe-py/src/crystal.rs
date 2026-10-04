//! Mechanical adapters for reciprocal geometry and exact reflection symmetry.

use crate::bindings::PyStructure;
use molframe::UnitCell;
use molframe::crystal::{AssemblyExt as _, AssemblyView};
use molframe::crystal::{
    CellTransform, ReflectionSymmetry, SymmetryExt as _, SymmetryOperation, SymmetrySet,
    space_group_by_hall, space_group_by_hermann_mauguin, space_group_setting, space_group_settings,
};
use numpy::{Complex64, PyArray1, PyArray2, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
use pyo3::prelude::*;

#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "UnitCell",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PyUnitCell(pub(crate) CellTransform, UnitCell);

impl PyUnitCell {
    /// The cell parameters.
    pub(crate) const fn parameters(&self) -> UnitCell {
        self.1
    }

    /// The cell a structure carries, as its conversion and its parameters.
    pub(crate) fn from_cell(cell: &UnitCell) -> PyResult<Self> {
        CellTransform::new(cell)
            .map(|transform| Self(transform, *cell))
            .map_err(crate::error::kernel)
    }
}

#[pymethods]
impl PyUnitCell {
    #[new]
    fn new(lengths: [f64; 3], angles: [f64; 3]) -> PyResult<Self> {
        Self::from_cell(&UnitCell { lengths, angles })
    }

    /// Edge lengths in ångström.
    #[getter]
    const fn lengths(&self) -> [f64; 3] {
        self.1.lengths
    }

    /// Angles in degrees.
    #[getter]
    const fn angles(&self) -> [f64; 3] {
        self.1.angles
    }

    #[pyo3(name = "to_cartesian")]
    fn fractional_to_cartesian(&self, fractional: [f64; 3]) -> [f64; 3] {
        self.0.to_cartesian(fractional)
    }

    #[pyo3(name = "to_fractional")]
    fn cartesian_to_fractional(&self, cartesian: [f64; 3]) -> [f64; 3] {
        self.0.to_fractional(cartesian)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.1 == other.1
    }

    fn __repr__(&self) -> String {
        format!(
            "UnitCell(lengths={:?}, angles={:?})",
            self.1.lengths, self.1.angles
        )
    }

    fn reciprocal_vector(&self, hkl: [i32; 3]) -> [f64; 3] {
        self.0.reciprocal_vector(hkl)
    }

    fn reciprocal_spacing_squared(&self, hkl: [i32; 3]) -> f64 {
        self.0.reciprocal_spacing_squared(hkl)
    }

    fn d_spacing(&self, hkl: [i32; 3]) -> f64 {
        self.0.d_spacing(hkl)
    }
}

#[derive(Clone, Debug)]
#[pyclass(
    name = "SpaceGroup",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PySpaceGroup(pub(crate) SymmetrySet);

#[pymethods]
impl PySpaceGroup {
    #[new]
    fn new(hall_number: u16) -> PyResult<Self> {
        space_group_setting(hall_number)
            .map(|setting| Self(setting.symmetry_set()))
            .map_err(crate::error::kernel)
    }

    /// The group a Hermann–Mauguin symbol names, short or full; case, spaces and screw-axis
    /// underscores do not matter, so a `CRYST1` spelling (`P 21 21 21`) resolves. A symbol of
    /// several settings gives the standard one.
    #[staticmethod]
    fn from_hermann_mauguin(symbol: &str) -> PyResult<Self> {
        space_group_by_hermann_mauguin(symbol)
            .map(|setting| Self(setting.symmetry_set()))
            .map_err(crate::error::kernel)
    }

    /// The setting a Hall symbol names.
    #[staticmethod]
    fn from_hall(symbol: &str) -> PyResult<Self> {
        space_group_by_hall(symbol)
            .map(|setting| Self(setting.symmetry_set()))
            .map_err(crate::error::kernel)
    }

    /// Every Hall setting of one International Tables type (`1..=230`).
    #[staticmethod]
    fn settings(international_number: u16) -> PyResult<Vec<Self>> {
        space_group_settings(international_number)
            .map(|found| {
                found
                    .iter()
                    .map(|setting| Self(setting.symmetry_set()))
                    .collect()
            })
            .map_err(crate::error::kernel)
    }

    /// The International Tables number, `1..=230`.
    #[getter]
    fn international_number(&self) -> Option<u16> {
        self.0.international_number
    }

    /// The Hermann–Mauguin symbol.
    #[getter]
    fn hermann_mauguin(&self) -> Option<&str> {
        self.0.hermann_mauguin.as_deref()
    }

    /// The unique-axis, origin or cell choice that distinguishes settings of one type.
    #[getter]
    fn choice(&self) -> Option<&str> {
        self.0.choice.as_deref()
    }

    /// The symmetry operations, each acting on fractional coordinates.
    #[getter]
    fn operations(&self) -> Vec<PySymmetryOperation> {
        self.0
            .operations()
            .iter()
            .cloned()
            .map(PySymmetryOperation)
            .collect()
    }

    fn __len__(&self) -> usize {
        self.0.operations().len()
    }

    fn __repr__(&self) -> String {
        match (&self.0.hermann_mauguin, self.0.international_number) {
            (Some(symbol), Some(number)) => format!("SpaceGroup({symbol}, number={number})"),
            _ => "SpaceGroup()".to_owned(),
        }
    }

    #[getter]
    fn hall_symbol(&self) -> Option<&str> {
        self.0.hall.as_deref()
    }

    #[getter]
    fn hall_number(&self) -> Option<u16> {
        self.0.hall_number
    }

    fn reflection_symmetry(&self, hkl: [i32; 3]) -> PyResult<PyReflectionSymmetry> {
        self.0
            .reflection_symmetry(hkl)
            .map(PyReflectionSymmetry)
            .map_err(crate::error::kernel)
    }
}

/// Three rows of three integers.
type RotationRows = ((i32, i32, i32), (i32, i32, i32), (i32, i32, i32));

/// One symmetry operation on fractional coordinates: `x' = W x + w`.
#[derive(Clone, Debug)]
#[pyclass(
    name = "SymmetryOperation",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PySymmetryOperation(SymmetryOperation);

#[pymethods]
impl PySymmetryOperation {
    /// The operation as written in the tables, such as `1/2-x,1/2+y,-z`.
    #[getter]
    fn expression(&self) -> String {
        self.0.to_string()
    }

    /// The integer rotation matrix `W`, row by row.
    #[getter]
    fn rotation(&self) -> RotationRows {
        let [a, b, c] = self.0.rotation;
        ((a[0], a[1], a[2]), (b[0], b[1], b[2]), (c[0], c[1], c[2]))
    }

    /// The exact fractional translation `w`.
    #[getter]
    fn translation(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.0.translation.map(molframe::crystal::Rational::as_f64);
        (x, y, z)
    }

    /// Whether this is the identity.
    #[getter]
    fn is_identity(&self) -> bool {
        self.0.is_identity()
    }

    /// Where the operation takes a fractional coordinate.
    fn apply(&self, fractional: [f64; 3]) -> [f64; 3] {
        self.0.apply_fractional(fractional)
    }

    fn __repr__(&self) -> String {
        format!("SymmetryOperation({})", self.0)
    }
}

/// The space group a structure carries, or `None` when it has none.
#[pyfunction]
fn space_group(py: Python<'_>, structure: &PyStructure) -> Option<PySpaceGroup> {
    let source = structure.inner.clone();
    py.detach(move || {
        source
            .engine()
            .symmetry_set()
            .map(|set| PySpaceGroup(set.clone()))
    })
}

#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "ReflectionSymmetry",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PyReflectionSymmetry(ReflectionSymmetry);

#[pymethods]
impl PyReflectionSymmetry {
    #[getter]
    fn centric(&self) -> bool {
        self.0.centric
    }
    #[getter]
    fn systematically_absent(&self) -> bool {
        self.0.systematically_absent
    }
    #[getter]
    fn epsilon_factor(&self) -> usize {
        self.0.epsilon_factor
    }
}

/// X-ray structure factors of the first model, one complex value per `(h, k, l)` row.
///
/// The cell, space group, occupancies and displacement parameters come from the
/// structure itself; the model is taken to be the asymmetric unit.
#[pyfunction]
fn structure_factors<'py>(
    py: Python<'py>,
    structure: &PyStructure,
    hkl: &Bound<'py, PyArray2<i32>>,
) -> PyResult<Bound<'py, PyArray1<Complex64>>> {
    let hkl = hkl.readonly();
    let shape = hkl.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(crate::error::value("hkl must have shape (n, 3)"));
    }
    let rows: Vec<[i32; 3]> = hkl
        .as_array()
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    let structure = structure.inner.clone();
    let values = py
        .detach(move || molframe::crystal::structure_factors(structure.engine(), &rows))
        .map_err(crate::error::kernel)?;
    let values: Vec<Complex64> = values
        .into_iter()
        .map(|value| Complex64::new(value.re, value.im))
        .collect();
    Ok(values.to_pyarray(py))
}

/// One placement of an assembly: a transform and the chains it applies to.
#[derive(Clone, Debug)]
#[pyclass(
    name = "AssemblyInstance",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
struct PyAssemblyInstance {
    matrix: [f64; 16],
    chains: Vec<String>,
}

#[pymethods]
impl PyAssemblyInstance {
    /// Column-major 4×4 affine matrix in ångström.
    #[getter]
    fn matrix(&self) -> Vec<f64> {
        self.matrix.to_vec()
    }

    /// `label_asym_id` of every chain the transform applies to.
    #[getter]
    fn chains(&self) -> Vec<String> {
        self.chains.clone()
    }

    fn __repr__(&self) -> String {
        format!("AssemblyInstance(chains={})", self.chains.join(","))
    }
}

/// Identifiers of the biological assemblies a structure declares.
#[pyfunction]
fn assemblies(structure: &PyStructure) -> Vec<String> {
    match structure.inner.engine().assembly_set() {
        Some(set) => set.assemblies().map(|each| each.id.to_string()).collect(),
        // A structure that declares no assemblies has none to list.
        None => Vec::new(),
    }
}

/// The placements that make up one biological assembly.
///
/// Each carries a transform and the chains it applies to; the transform of
/// every chain is the product of the operators the entry lists for it.
#[pyfunction]
fn assembly(
    py: Python<'_>,
    structure: &PyStructure,
    id: &str,
) -> PyResult<Vec<PyAssemblyInstance>> {
    let source = structure.inner.clone();
    py.detach(|| placements(&source, id))
}

fn placements(structure: &molframe::Structure, id: &str) -> PyResult<Vec<PyAssemblyInstance>> {
    let engine = structure.engine();
    let set = engine
        .assembly_set()
        .ok_or_else(|| crate::error::value("the structure declares no biological assemblies"))?;
    let view = AssemblyView::new(engine, set, id).map_err(crate::error::kernel)?;
    view.groups_by_transform()
        .into_iter()
        .map(|(transform, chains)| {
            let mut matrix = [0.0_f64; 16];
            for column in 0..3 {
                for row in 0..3 {
                    matrix[column * 4 + row] = transform.rotation[row][column];
                }
            }
            matrix[12..15].copy_from_slice(&transform.translation);
            matrix[15] = 1.0;
            let labels = chains
                .into_iter()
                .map(|chain| {
                    engine
                        .chain(chain)
                        .and_then(molframe::ChainRef::label)
                        .map(str::to_owned)
                        .ok_or_else(|| crate::error::value("an assembly chain has no label"))
                })
                .collect::<PyResult<Vec<_>>>()?;
            Ok(PyAssemblyInstance {
                matrix,
                chains: labels,
            })
        })
        .collect()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::crystal_links::register(module)?;
    crate::crystal_maps::register(module)?;
    crate::crystal_reflections::register(module)?;
    module.add_class::<PyAssemblyInstance>()?;
    module.add_function(wrap_pyfunction!(assemblies, module)?)?;
    module.add_function(wrap_pyfunction!(assembly, module)?)?;
    module.add_function(wrap_pyfunction!(structure_factors, module)?)?;
    module.add_class::<PyUnitCell>()?;
    module.add_class::<PySpaceGroup>()?;
    module.add_class::<PySymmetryOperation>()?;
    module.add_function(wrap_pyfunction!(space_group, module)?)?;
    module.add_class::<PyReflectionSymmetry>()?;
    crate::crystal_statistics::register(module)?;
    crate::crystal_reduction::register(module)
}
