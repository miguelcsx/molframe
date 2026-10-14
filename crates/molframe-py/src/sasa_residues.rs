//! Per-residue solvent-accessible surface area.
//!
//! Rust sums the per-atom areas of each residue under the policy's radii rule;
//! this layer moves the columns into read-only arrays.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::governed::run;
use crate::governed_field::sasa_options;
use crate::policy::{PyAnalysisPolicy, policy_of};
use numpy::{PyArray1, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyList;

/// Solvent-accessible area of every residue, one row per residue in chain order.
///
/// `chain`, `name` are lists of strings; `number`, `area` and `atoms` are
/// read-only arrays of the same length.
#[derive(Debug)]
#[pyclass(
    name = "ResidueAreas",
    frozen,
    skip_from_py_object,
    module = "molframe.analysis"
)]
struct PyResidueAreas {
    chain: Py<PyList>,
    name: Py<PyList>,
    number: Py<PyArray1<i32>>,
    area: Py<PyArray1<f64>>,
    atoms: Py<PyArray1<u64>>,
    excluded_hetero_residues: usize,
    skipped_atoms: usize,
    missing_radius: usize,
}

#[pymethods]
impl PyResidueAreas {
    /// Chain label of every residue.
    #[getter]
    fn chain(&self, py: Python<'_>) -> Py<PyList> {
        self.chain.clone_ref(py)
    }

    /// Residue name of every residue.
    #[getter]
    fn name(&self, py: Python<'_>) -> Py<PyList> {
        self.name.clone_ref(py)
    }

    /// Residue number of every residue.
    #[getter]
    fn number(&self, py: Python<'_>) -> Py<PyArray1<i32>> {
        self.number.clone_ref(py)
    }

    /// Accessible area of every residue in square angstroms: the sum of its atoms' areas.
    #[getter]
    fn area(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.area.clone_ref(py)
    }

    /// Measured atoms of every residue.
    #[getter]
    fn atoms(&self, py: Python<'_>) -> Py<PyArray1<u64>> {
        self.atoms.clone_ref(py)
    }

    /// Heteroatom residues the options left out.
    #[getter]
    fn excluded_hetero_residues(&self) -> usize {
        self.excluded_hetero_residues
    }

    /// Atoms left out: hydrogens, other conformers and heteroatoms not selected.
    #[getter]
    fn skipped_atoms(&self) -> usize {
        self.skipped_atoms
    }

    /// Included atoms the radius table has no radius for.
    #[getter]
    fn missing_radius(&self) -> usize {
        self.missing_radius
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.area.bind(py).len()
    }
}

fn widen(count: usize) -> u64 {
    match u64::try_from(count) {
        Ok(count) => count,
        Err(_) => u64::MAX,
    }
}

fn readonly<T: numpy::Element>(py: Python<'_>, values: Vec<T>) -> Py<PyArray1<T>> {
    let array = PyArray1::from_vec(py, values);
    array.readwrite().make_nonwriteable();
    array.unbind()
}

fn convert(py: Python<'_>, areas: molframe::analysis::ResidueAreas) -> PyResult<Py<PyAny>> {
    let rows = areas.residues.len();
    let mut chain = Vec::with_capacity(rows);
    let mut name = Vec::with_capacity(rows);
    let mut number = Vec::with_capacity(rows);
    let mut area = Vec::with_capacity(rows);
    let mut atoms = Vec::with_capacity(rows);
    for residue in areas.residues {
        chain.push(residue.chain);
        name.push(residue.name);
        number.push(residue.number);
        area.push(residue.area);
        atoms.push(widen(residue.atoms));
    }
    let table = PyResidueAreas {
        chain: PyList::new(py, chain)?.unbind(),
        name: PyList::new(py, name)?.unbind(),
        number: readonly(py, number),
        area: readonly(py, area),
        atoms: readonly(py, atoms),
        excluded_hetero_residues: areas.excluded_hetero_residues,
        skipped_atoms: areas.skipped_atoms,
        missing_radius: areas.missing_radius,
    };
    Ok(Py::new(py, table)?.into_any())
}

/// Solvent-accessible surface area per residue, as the sum of its atoms' areas.
///
/// Options and radii follow `sasa`; the result records the radius table that applied.
#[pyfunction]
#[pyo3(signature = (structure, *, probe=None, method=None, n_points=None, hetero=None, metals=None, radii=None, policy=None, context=None))]
pub(crate) fn residue_sasa(
    py: Python<'_>,
    structure: &PyStructure,
    probe: Option<f32>,
    method: Option<&str>,
    n_points: Option<u16>,
    hetero: Option<&str>,
    metals: Option<&str>,
    radii: Option<&str>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let (selection, source) = crate::tables::radii_under(radii, policy.as_deref())?;
    let options = sasa_options(method, n_points, probe, hetero, metals, &source)?;
    let kernel = molframe::analysis::residue_sasa_kernel(options, selection);
    run(py, structure, &policy_of(policy), &kernel, convert, context)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyResidueAreas>()?;
    module.add_function(wrap_pyfunction!(residue_sasa, module)?)
}
