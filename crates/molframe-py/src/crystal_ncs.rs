//! Non-crystallographic symmetry: operators, the chains they relate, and the
//! structure with the missing copies placed. The work happens in the kernel;
//! this module only hands its results across.

use crate::bindings::PyStructure;
use molframe::crystal::{
    DEFAULT_NCS_PAIR_LIMIT, NcsBindOptions, NcsChainPair, NcsCode, NcsExt as _, NcsOperator,
};
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;

/// One NCS operator.
#[pyclass(name = "NcsOperator", module = "molframe.crystal", skip_from_py_object)]
struct PyNcsOperator {
    id: String,
    code: &'static str,
    details: Option<String>,
    rotation: [f64; 9],
    translation: [f64; 3],
}

#[pymethods]
impl PyNcsOperator {
    /// Identifier in the source file.
    #[getter]
    fn id(&self) -> String {
        self.id.clone()
    }

    /// `given` when the copy is in the coordinates, `generate` when it is not.
    #[getter]
    const fn code(&self) -> &'static str {
        self.code
    }

    /// Free-text description, when the source gave one.
    #[getter]
    fn details(&self) -> Option<String> {
        self.details.clone()
    }

    /// Row-major 3x3 matrix.
    #[getter]
    fn rotation<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        PyArray1::from_slice(py, &self.rotation).reshape([3, 3])
    }

    /// Translation in ångström, shape (3,).
    #[getter]
    fn translation<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &self.translation)
    }

    fn __repr__(&self) -> String {
        format!("NcsOperator(id={:?}, code={:?})", self.id, self.code)
    }
}

/// One operator bound to the chain it maps and the chain it produces.
#[pyclass(
    name = "NcsChainPair",
    module = "molframe.crystal",
    skip_from_py_object
)]
struct PyNcsChainPair {
    operator: String,
    source: String,
    target: String,
    rmsd: Option<f64>,
    atoms: usize,
    declared: bool,
}

#[pymethods]
impl PyNcsChainPair {
    /// Operator identifier.
    #[getter]
    fn operator(&self) -> String {
        self.operator.clone()
    }

    /// `label_asym_id` the operator is applied to.
    #[getter]
    fn source(&self) -> String {
        self.source.clone()
    }

    /// `label_asym_id` the operator produces.
    #[getter]
    fn target(&self) -> String {
        self.target.clone()
    }

    /// Backbone RMSD of the mapped source against the target, in ångström.
    #[getter]
    const fn rmsd(&self) -> Option<f64> {
        self.rmsd
    }

    /// Number of atom pairs the RMSD covers.
    #[getter]
    const fn atoms(&self) -> usize {
        self.atoms
    }

    /// Whether the file declared the pairing rather than geometry finding it.
    #[getter]
    const fn declared(&self) -> bool {
        self.declared
    }

    fn __repr__(&self) -> String {
        format!(
            "NcsChainPair(operator={:?}, source={:?}, target={:?})",
            self.operator, self.source, self.target
        )
    }
}

fn operator_of(operator: &NcsOperator) -> PyNcsOperator {
    let matrix = operator.transform.matrix;
    let mut rotation = [0.0; 9];
    for (row, values) in matrix.iter().enumerate() {
        rotation[row * 3..row * 3 + 3].copy_from_slice(values);
    }
    PyNcsOperator {
        id: operator.id.to_string(),
        code: match operator.code {
            NcsCode::Given => "given",
            NcsCode::Generate => "generate",
        },
        details: operator.details.as_deref().map(str::to_owned),
        rotation,
        translation: operator.transform.translation,
    }
}

/// The NCS operators the structure records; empty when it records none.
#[pyfunction]
fn ncs_operators(structure: &PyStructure) -> Vec<PyNcsOperator> {
    structure
        .inner
        .ncs_set()
        .map(|set| set.operators().map(operator_of).collect())
        .unwrap_or_default()
}

fn label_of(structure: &molframe::Structure, chain: molframe::ChainIndex) -> PyResult<String> {
    structure
        .engine()
        .chain(chain)
        .and_then(molframe::ChainRef::label)
        .map(str::to_owned)
        .ok_or_else(|| crate::error::value("an NCS chain has no label"))
}

/// The chain pairs the NCS operators relate.
///
/// A file that declares domain pairings binds exactly those; otherwise each
/// `given` operator is tested against chains of identical sequence and kept
/// where the mapped backbone lies within `max_rmsd` of the other chain.
#[pyfunction]
#[pyo3(signature = (structure, max_rmsd = 2.0, include_generate = false, pair_limit = DEFAULT_NCS_PAIR_LIMIT))]
fn ncs_chain_pairs(
    py: Python<'_>,
    structure: &PyStructure,
    max_rmsd: f64,
    include_generate: bool,
    pair_limit: usize,
) -> PyResult<Vec<PyNcsChainPair>> {
    let options = NcsBindOptions {
        max_rmsd,
        include_generate,
        pair_limit,
    };
    let source = structure.inner.clone();
    let pairs = py
        .detach(|| source.ncs_chain_pairs(&options))
        .map_err(crate::error::kernel)?;
    pairs
        .iter()
        .map(|pair: &NcsChainPair| {
            Ok(PyNcsChainPair {
                operator: pair.operator.to_string(),
                source: label_of(&source, pair.source)?,
                target: label_of(&source, pair.target)?,
                rmsd: pair.rmsd,
                atoms: pair.atoms,
                declared: pair.declared,
            })
        })
        .collect()
}

/// The structure with every `generate` operator applied to all of its chains.
///
/// The deposited coordinates come first, followed by one transformed copy of
/// the model per operator; a structure without such operators is returned
/// unchanged.
#[pyfunction]
fn ncs_structure(py: Python<'_>, structure: &PyStructure) -> PyResult<PyStructure> {
    let source = structure.inner.clone();
    let built = py
        .detach(|| source.ncs_expanded())
        .map_err(|findings| crate::error::from_findings(&findings, "NCS expansion failed"))?;
    Ok(PyStructure::new(molframe::Structure::from(built)))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNcsOperator>()?;
    module.add_class::<PyNcsChainPair>()?;
    module.add_function(wrap_pyfunction!(ncs_operators, module)?)?;
    module.add_function(wrap_pyfunction!(ncs_chain_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(ncs_structure, module)?)
}
