//! Residue-aware lDDT of a model structure against a reference structure.

use crate::bindings::PyStructure;
use crate::execution::PyExecutionContext;
use molframe::compare::MappedCompareError;
use molframe::compare::{
    EmptyLddtPolicy, ResidueLddt, ResidueLddtOptions, SymmetryPolicy, pair_residues_by_identifier,
    residue_lddt,
};
use numpy::PyArray1;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;

/// A residue-aware lDDT and the choices it was computed under.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Lddt",
    frozen,
    skip_from_py_object,
    module = "molframe.compare"
)]
pub(crate) struct PyLddt {
    inner: ResidueLddt,
    options: ResidueLddtOptions,
}

#[pymethods]
impl PyLddt {
    /// Fraction of tolerance tests passed over every tested distance, in `[0, 1]`.
    #[getter]
    const fn score(&self) -> f64 {
        self.inner.score
    }

    /// Reference distances tested, each pair counted once.
    #[getter]
    const fn considered(&self) -> u64 {
        self.inner.considered
    }

    /// Residues whose interchangeable atoms were renamed to a better naming.
    #[getter]
    const fn renamed_residues(&self) -> usize {
        self.inner.renamed_residues
    }

    /// Residues whose symmetry exceeded `automorphism_limit` and was left alone.
    #[getter]
    const fn unresolved_symmetry(&self) -> usize {
        self.inner.unresolved_symmetry
    }

    /// The score of each reference residue, in reference order; `nan` where
    /// the residue has no tested distance.
    #[getter]
    fn residue_scores<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let scores = self
            .inner
            .residues
            .iter()
            .map(|residue| match residue.score {
                Some(score) => score,
                None => f64::NAN,
            })
            .collect();
        PyArray1::from_vec(py, scores)
    }

    /// Reference residue indices, in the order of `residue_scores`.
    #[getter]
    fn reference_residues<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        let rows = self.inner.residues.iter().map(|r| r.reference.get());
        PyArray1::from_vec(py, rows.collect())
    }

    /// The model residue each reference residue was compared with, or `-1`
    /// where the model has none; such a residue's distances count as lost.
    #[getter]
    fn model_residues<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        let rows = self
            .inner
            .residues
            .iter()
            .map(|residue| match residue.model {
                Some(model) => i64::from(model.get()),
                None => -1,
            });
        PyArray1::from_vec(py, rows.collect())
    }

    /// The parameters and fixed conventions the score was computed under.
    #[getter]
    fn parameters<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let options = &self.options;
        let parameters = PyDict::new(py);
        parameters.set_item("inclusion_radius", options.inclusion_radius)?;
        parameters.set_item("tolerances", options.tolerances.to_vec())?;
        parameters.set_item("sequence_separation", options.sequence_separation)?;
        let symmetry = match options.symmetry {
            SymmetryPolicy::Disabled => "disabled",
            SymmetryPolicy::Resolve { .. } => "topology-only",
        };
        parameters.set_item("symmetry", symmetry)?;
        parameters.set_item("reference_pairs", "heavy-atom-inter-residue")?;
        parameters.set_item("missing_model_atoms", "not-preserved")?;
        parameters.set_item("pair_counting", "both-sides")?;
        Ok(parameters)
    }

    fn __repr__(&self) -> String {
        format!(
            "Lddt(score={:.4}, considered={}, renamed_residues={})",
            self.inner.score, self.inner.considered, self.inner.renamed_residues
        )
    }
}

/// lDDT of `model` against `reference`, residue by residue.
///
/// Reference distances up to `inclusion_radius` Å between atoms of different
/// residues are tested; each counts as preserved at every tolerance it stays
/// within, and the score is the preserved fraction. A reference atom the model
/// lacks, or a whole reference residue the model lacks, counts as not
/// preserved. Residues correspond by chain and sequence position under
/// `chain_names`, and only when they have the same component name.
/// `components` is a Chemical Component Dictionary file; it says which atoms
/// of a residue are interchangeable. With `symmetry` they are renamed to the
/// naming that preserves the most reference distances, so a model that calls a
/// carboxylate oxygen the other name is not penalised. Residues more than
/// `sequence_separation` apart in one chain are tested only past that distance.
#[pyfunction]
#[pyo3(signature = (
    model,
    reference,
    *,
    components,
    components_version,
    inclusion_radius=15.0,
    tolerances=vec![0.5, 1.0, 2.0, 4.0],
    sequence_separation=0,
    symmetry=true,
    automorphism_limit=64,
    chain_names="label",
    context=None,
))]
pub(crate) fn lddt(
    py: Python<'_>,
    model: &PyStructure,
    reference: &PyStructure,
    components: PathBuf,
    components_version: &str,
    inclusion_radius: f64,
    tolerances: Vec<f64>,
    sequence_separation: u32,
    symmetry: bool,
    automorphism_limit: usize,
    chain_names: &str,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyLddt> {
    let namespace = chain_names.parse().map_err(crate::error::kernel)?;
    let (provider, _) = molframe::read_component_dictionary(
        components,
        molframe::DictionaryVersion::new(components_version),
    )
    .map_err(|findings| crate::bindings::findings_error(&findings))?;
    let options = ResidueLddtOptions {
        inclusion_radius,
        tolerances: tolerances.into_boxed_slice(),
        sequence_separation,
        symmetry: if symmetry {
            SymmetryPolicy::Resolve { automorphism_limit }
        } else {
            SymmetryPolicy::Disabled
        },
        empty_policy: EmptyLddtPolicy::Perfect,
    };
    let (model, reference) = (model.inner.clone(), reference.inner.clone());
    let found = crate::execution::run(py, context, |context| {
        let pairs =
            pair_residues_by_identifier(reference.engine(), model.engine(), &provider, namespace)
                .map_err(MappedCompareError::Mapping)?;
        residue_lddt(
            reference.engine(),
            model.engine(),
            &pairs,
            &provider,
            &options,
            context,
        )
    })?
    .map_err(crate::error::kernel)?;
    Ok(PyLddt {
        inner: found,
        options,
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyLddt>()?;
    module.add_function(wrap_pyfunction!(lddt, module)?)
}
