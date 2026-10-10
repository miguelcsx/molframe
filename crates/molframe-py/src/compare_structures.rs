//! Comparing structures and models: local-distance scores, docking quality
//! and quaternary structure, with and without a chain and atom mapping.

use crate::bindings::PyStructure;
use crate::execution::PyExecutionContext;
use crate::sequence::PyScoring;
use molframe::compare::{
    DockQ, DockQOptions, EmptyLddtPolicy, LddtOptions, MappingOptions, QsOptions,
};
use molframe::{CompareExt, Namespace};
use numpy::{PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use std::path::PathBuf;

fn namespace(name: &str) -> PyResult<Namespace> {
    name.parse().map_err(crate::error::kernel)
}

/// The local distance difference test of a model against a reference.
///
/// Every reference distance up to `inclusion_radius` Å is tested, and a model
/// distance counts as preserved at each tolerance in `tolerances` that it stays
/// within; the score is the preserved fraction averaged over the tolerances.
/// A reference with no such distance scores `1.0`: there is nothing to get
/// wrong. The defaults are the standard 15 Å radius and 0.5/1/2/4 Å tolerances.
#[pyfunction]
#[pyo3(signature = (
    model,
    reference,
    *,
    inclusion_radius=15.0,
    tolerances=vec![0.5, 1.0, 2.0, 4.0],
    minimum_reference_distance=0.0,
    context=None,
))]
pub(crate) fn lddt(
    py: Python<'_>,
    model: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
    inclusion_radius: f64,
    tolerances: Vec<f64>,
    minimum_reference_distance: f64,
    context: Option<&PyExecutionContext>,
) -> PyResult<f64> {
    let (model, reference) = (model.readonly(), reference.readonly());
    let (model, reference) = (
        crate::bindings::coordinates(&model)?,
        crate::bindings::coordinates(&reference)?,
    );
    let options = LddtOptions {
        inclusion_radius,
        minimum_reference_distance,
        tolerances: tolerances.into_boxed_slice(),
        empty_policy: EmptyLddtPolicy::Perfect,
    };
    crate::execution::run(py, context, |context| {
        molframe::compare::lddt_with_options(model, reference, &options, context)
    })?
    .map_err(crate::error::kernel)
}

/// A docking score and its three components.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "DockQ",
    frozen,
    skip_from_py_object,
    module = "molframe.compare"
)]
pub(crate) struct PyDockQ {
    inner: DockQ,
}

#[pymethods]
impl PyDockQ {
    /// Fraction of native interface contacts the model keeps.
    #[getter]
    const fn fnat(&self) -> f64 {
        self.inner.fnat
    }

    /// Ligand RMSD after superposing the receptor, in ångström.
    #[getter]
    const fn ligand_rmsd(&self) -> f64 {
        self.inner.ligand_rmsd
    }

    /// Interface RMSD, in ångström.
    #[getter]
    const fn interface_rmsd(&self) -> f64 {
        self.inner.interface_rmsd
    }

    /// The combined score in `(0, 1]`.
    #[getter]
    const fn score(&self) -> f64 {
        self.inner.score
    }

    fn __repr__(&self) -> String {
        format!(
            "DockQ(score={:.4}, fnat={:.4}, ligand_rmsd={:.3}, interface_rmsd={:.3})",
            self.inner.score, self.inner.fnat, self.inner.ligand_rmsd, self.inner.interface_rmsd
        )
    }
}

/// What corresponded when chains and atoms were mapped before scoring.
#[derive(Clone, Debug)]
#[pyclass(
    name = "MappedComparison",
    frozen,
    skip_from_py_object,
    module = "molframe.compare"
)]
pub(crate) struct PyMappedComparison {
    chains: Vec<(String, String, f64)>,
    atoms: usize,
    swapped_residues: usize,
    status: &'static str,
    chain_assignments_tried: usize,
}

#[pymethods]
impl PyMappedComparison {
    /// `(native chain, model chain, sequence identity)` for the scored chains.
    #[getter]
    fn chains(&self) -> Vec<(String, String, f64)> {
        self.chains.clone()
    }

    /// How many atom pairs were compared.
    #[getter]
    const fn atoms(&self) -> usize {
        self.atoms
    }

    /// Residue pairs whose chemically equivalent atoms were swapped to fit better.
    #[getter]
    const fn swapped_residues(&self) -> usize {
        self.swapped_residues
    }

    /// `"complete"`, or `"ambiguous"` when tied atom correspondences exist.
    #[getter]
    const fn status(&self) -> &'static str {
        self.status
    }

    /// How many chain assignments of equal sequence identity were scored; more
    /// than one means the sequences alone could not say which model chain
    /// answers to which native chain, and the best-scoring assignment is reported.
    #[getter]
    const fn chain_assignments_tried(&self) -> usize {
        self.chain_assignments_tried
    }
}

impl PyMappedComparison {
    fn from_native(found: &molframe::compare::MappedComparison) -> Self {
        Self {
            chains: found
                .chains
                .iter()
                .map(|chain| {
                    (
                        chain.reference.clone(),
                        chain.target.clone(),
                        chain.identity,
                    )
                })
                .collect(),
            atoms: found.atoms.matches().len(),
            swapped_residues: found.swapped_residues,
            status: match found.status {
                molframe::Status::Complete => "complete",
                molframe::Status::Partial => "partial",
                molframe::Status::Ambiguous => "ambiguous",
                molframe::Status::Indeterminate => "indeterminate",
                _ => "unknown",
            },
            chain_assignments_tried: found.chain_assignments_tried,
        }
    }
}

/// `DockQ` of a docking `model` against the `native` complex.
///
/// The structures must have the same atoms in the same order; use
/// `mapped_dockq` when chains are renamed or atoms ordered differently. The
/// defaults are the published `DockQ` constants: a 5 Å contact distance and 8.5 Å
/// and 1.5 Å length scales for the ligand and interface RMSD.
#[pyfunction]
#[pyo3(signature = (
    model,
    native,
    *,
    receptor,
    ligand,
    contact_distance=5.0,
    ligand_scale=8.5,
    interface_scale=1.5,
    chain_names="label",
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn dockq(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    receptor: &str,
    ligand: &str,
    contact_distance: f32,
    ligand_scale: f64,
    interface_scale: f64,
    chain_names: &str,
) -> PyResult<PyDockQ> {
    let namespace = namespace(chain_names)?;
    let options = DockQOptions {
        contact_distance,
        ligand_scale,
        interface_scale,
        ..DockQOptions::published()
    };
    let (model, native) = (model.inner.clone(), native.inner.clone());
    py.detach(|| model.dockq_in_namespace(&native, receptor, ligand, namespace, options))
        .map(|inner| PyDockQ { inner })
        .map_err(crate::error::kernel)
}

/// Residue-level, distance-weighted QS score between two chains (CB, or CA for glycine; a
/// `contact_distance` of 12 reproduces the published cutoff).
#[pyfunction]
#[pyo3(signature = (
    model,
    native,
    *,
    first_chain,
    second_chain,
    contact_distance=5.0,
    chain_names="label",
))]
pub(crate) fn qs_score(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    first_chain: &str,
    second_chain: &str,
    contact_distance: f32,
    chain_names: &str,
) -> PyResult<f64> {
    let namespace = namespace(chain_names)?;
    let options = QsOptions::standard(contact_distance);
    let (model, native) = (model.inner.clone(), native.inner.clone());
    py.detach(|| {
        model.qs_score_in_namespace(&native, first_chain, second_chain, namespace, options)
    })
    .map_err(crate::error::kernel)
}

/// The mapping decisions a call states: component definitions and alignment.
struct Mapping {
    provider: Box<dyn molframe::chemistry::ComponentProvider + Send + Sync>,
    namespace: Namespace,
    scoring: molframe::sequence::Scoring,
    min_identity: f64,
    automorphism_limit: usize,
}

impl Mapping {
    fn load(
        components: &std::path::Path,
        components_version: &str,
        chain_names: &str,
        scoring: &PyScoring,
        min_identity: f64,
        automorphism_limit: usize,
    ) -> PyResult<Self> {
        let (provider, _) = molframe::read_component_dictionary(
            components,
            molframe::DictionaryVersion::new(components_version),
        )
        .map_err(|findings| crate::bindings::findings_error(&findings))?;
        Ok(Self {
            provider: Box::new(provider),
            namespace: namespace(chain_names)?,
            scoring: scoring.0,
            min_identity,
            automorphism_limit,
        })
    }

    fn options(&self) -> MappingOptions<'_> {
        MappingOptions {
            provider: self.provider.as_ref(),
            namespace: self.namespace,
            scoring: self.scoring,
            min_identity: self.min_identity,
            automorphism_limit: self.automorphism_limit,
        }
    }
}

fn compare_error(error: &molframe::compare::MappedCompareError) -> PyErr {
    crate::error::kernel(error)
}

/// `DockQ` after matching chains by sequence and atoms by residue alignment.
///
/// `components` is a Chemical Component Dictionary file giving each residue's
/// sequence code and chemically equivalent atoms; `scoring` and `min_identity`
/// state how chains correspond. Chains of equal identity (the halves of a
/// homodimer) are resolved by scoring every assignment of the named pair, and
/// `comparison.chain_assignments_tried` says how many there were.
#[pyfunction]
#[pyo3(signature = (
    model,
    native,
    *,
    receptor,
    ligand,
    components,
    components_version,
    scoring,
    min_identity,
    automorphism_limit,
    contact_distance=5.0,
    ligand_scale=8.5,
    interface_scale=1.5,
    chain_names="label",
))]
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
pub(crate) fn mapped_dockq(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    receptor: &str,
    ligand: &str,
    components: PathBuf,
    components_version: &str,
    scoring: &PyScoring,
    min_identity: f64,
    automorphism_limit: usize,
    contact_distance: f32,
    ligand_scale: f64,
    interface_scale: f64,
    chain_names: &str,
) -> PyResult<(PyDockQ, PyMappedComparison)> {
    let mapping = Mapping::load(
        &components,
        components_version,
        chain_names,
        scoring,
        min_identity,
        automorphism_limit,
    )?;
    let options = DockQOptions {
        contact_distance,
        ligand_scale,
        interface_scale,
        ..DockQOptions::published()
    };
    let (model, native) = (model.inner.clone(), native.inner.clone());
    py.detach(|| model.mapped_dockq(&native, receptor, ligand, &mapping.options(), options))
        .map(|(score, found)| {
            (
                PyDockQ { inner: score },
                PyMappedComparison::from_native(&found),
            )
        })
        .map_err(|error| compare_error(&error))
}

/// QS score after matching chains by sequence and atoms by residue alignment.
#[pyfunction]
#[pyo3(signature = (
    model,
    native,
    *,
    first_chain,
    second_chain,
    components,
    components_version,
    scoring,
    min_identity,
    automorphism_limit,
    contact_distance=5.0,
    chain_names="label",
))]
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
pub(crate) fn mapped_qs_score(
    py: Python<'_>,
    model: &PyStructure,
    native: &PyStructure,
    first_chain: &str,
    second_chain: &str,
    components: PathBuf,
    components_version: &str,
    scoring: &PyScoring,
    min_identity: f64,
    automorphism_limit: usize,
    contact_distance: f32,
    chain_names: &str,
) -> PyResult<(f64, PyMappedComparison)> {
    let mapping = Mapping::load(
        &components,
        components_version,
        chain_names,
        scoring,
        min_identity,
        automorphism_limit,
    )?;
    let options = QsOptions::standard(contact_distance);
    let (model, native) = (model.inner.clone(), native.inner.clone());
    py.detach(|| {
        model.mapped_qs_score(
            &native,
            first_chain,
            second_chain,
            &mapping.options(),
            options,
        )
    })
    .map(|(score, found)| (score, PyMappedComparison::from_native(&found)))
    .map_err(|error| compare_error(&error))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDockQ>()?;
    module.add_class::<PyMappedComparison>()?;
    module.add_function(wrap_pyfunction!(lddt, module)?)?;
    module.add_function(wrap_pyfunction!(dockq, module)?)?;
    module.add_function(wrap_pyfunction!(qs_score, module)?)?;
    module.add_function(wrap_pyfunction!(mapped_dockq, module)?)?;
    module.add_function(wrap_pyfunction!(mapped_qs_score, module)?)
}
