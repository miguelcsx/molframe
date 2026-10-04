//! Governed validation checks: bond geometry, peptide planarity, occupancy,
//! completeness and displacement factors.
//!
//! Every tolerance is a required keyword. Each result carries status, coverage
//! and provenance, and a check that assessed less than it intended says so in
//! its coverage rather than reporting a clean pass.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::execution::PyExecutionContext;
use crate::governed::run;
use crate::policy::{PyAnalysisPolicy, policy_of};
use crate::table::TableBuilder;
use molframe::geometry::EigenOptions;
use molframe::validation::{AltlocOccupancyOptions, PlanarityOptions, QualityIssue};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

fn atoms(column: impl Iterator<Item = molframe::AtomIndex>) -> Vec<u32> {
    column.map(molframe::AtomIndex::get).collect()
}

fn residues(column: impl Iterator<Item = molframe::ResidueIndex>) -> Vec<u32> {
    column.map(molframe::ResidueIndex::get).collect()
}

/// Bonds whose length differs from the covalent expectation by more than
/// `tolerance` ångström: `atom_a`, `atom_b`, `observed`, `expected`, `deviation`.
#[pyfunction]
#[pyo3(signature = (structure, *, tolerance, policy=None, context=None))]
pub(crate) fn bond_length_deviations(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::bond_length_deviations_kernel(tolerance);
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let observed: Vec<f32> = rows.iter().map(|row| row.observed).collect();
            let expected: Vec<f32> = rows.iter().map(|row| row.expected).collect();
            let deviation: Vec<f32> = rows.iter().map(|row| row.deviation).collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("atom_a", &atoms(rows.iter().map(|row| row.atom_a)))
                .indices("atom_b", &atoms(rows.iter().map(|row| row.atom_b)))
                .single("observed", &observed)
                .single("expected", &expected)
                .single("deviation", &deviation)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Peptide bonds whose ω torsion is within `threshold_degrees` of zero:
/// `residue` (the first of the pair) and `omega` in degrees.
#[pyfunction]
#[pyo3(signature = (structure, *, threshold_degrees, policy=None, context=None))]
pub(crate) fn cis_peptides(
    py: Python<'_>,
    structure: &PyStructure,
    threshold_degrees: f64,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::cis_peptides_kernel(threshold_degrees);
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let omega: Vec<f64> = rows.iter().map(|row| row.omega).collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("residue", &residues(rows.iter().map(|row| row.residue)))
                .double("omega", &omega)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Atoms with a zero or out-of-range occupancy or a negative B factor.
///
/// `issue` is `0` for zero occupancy, `1` for an occupancy outside `[0, 1]`
/// and `2` for a negative B factor.
#[pyfunction]
#[pyo3(signature = (structure, *, policy=None, context=None))]
pub(crate) fn quality_flags(
    py: Python<'_>,
    structure: &PyStructure,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::quality_flags_kernel();
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let issue: Vec<u32> = rows
                .iter()
                .map(|row| match row.issue {
                    QualityIssue::ZeroOccupancy => 0,
                    QualityIssue::OccupancyOutOfRange => 1,
                    QualityIssue::NegativeBFactor => 2,
                })
                .collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("atom", &atoms(rows.iter().map(|row| row.atom)))
                .indices("issue", &issue)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Atoms with more bonds than their element allows: `atom`, `bonds`, `maximum`.
#[pyfunction]
#[pyo3(signature = (structure, *, policy=None, context=None))]
pub(crate) fn valence(
    py: Python<'_>,
    structure: &PyStructure,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::valence_kernel();
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let count = |pick: fn(&molframe::validation::ValenceError) -> usize| -> Vec<u32> {
                rows.iter()
                    .map(|row| match u32::try_from(pick(row)) {
                        Ok(value) => value,
                        Err(_) => u32::MAX,
                    })
                    .collect()
            };
            let table = TableBuilder::new(py, rows.len())
                .indices("atom", &atoms(rows.iter().map(|row| row.atom)))
                .indices("bonds", &count(|row| row.bonds))
                .indices("maximum", &count(|row| row.maximum))
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Residues of aromatic rings that deviate from their plane by more than
/// `max_deviation` ångström: `residue` and `deviation`.
#[pyfunction]
#[pyo3(signature = (structure, *, max_deviation, policy=None, context=None))]
pub(crate) fn planarity(
    py: Python<'_>,
    structure: &PyStructure,
    max_deviation: f64,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::planarity_kernel(PlanarityOptions {
        maximum_deviation: max_deviation,
        plane_fit: EigenOptions::standard(),
    });
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, rows| {
            let deviation: Vec<f64> = rows.iter().map(|row| row.deviation).collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("residue", &residues(rows.iter().map(|row| row.residue)))
                .double("deviation", &deviation)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Per polymer chain, how many residues were modelled against the entity's
/// canonical sequence and which canonical positions are missing.
///
/// The value is a list of dictionaries with `chain`, `observed`, `canonical`
/// and `missing` (a list of `(canonical position, component)` pairs). Chain
/// names follow the policy's identifier namespace (`label` or `auth`).
#[pyfunction]
#[pyo3(signature = (structure, *, policy=None, context=None))]
pub(crate) fn completeness(
    py: Python<'_>,
    structure: &PyStructure,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::completeness_kernel();
    let source = structure.inner.clone();
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, chains| {
            let list = PyList::empty(py);
            for chain in chains {
                let entry = PyDict::new(py);
                entry.set_item("chain", chain.chain)?;
                entry.set_item("observed", chain.observed)?;
                entry.set_item("canonical", chain.canonical)?;
                let missing: Vec<(u32, String)> = chain
                    .missing
                    .iter()
                    .map(|residue| {
                        let name = match source.engine().resolve(residue.component) {
                            Some(name) => name.to_owned(),
                            None => String::new(),
                        };
                        (residue.canonical_position, name)
                    })
                    .collect();
                entry.set_item("missing", missing)?;
                list.append(entry)?;
            }
            Ok(list.into_any().unbind())
        },
        context,
    )
}

/// Alternate-location groups whose occupancies do not sum to `expected_sum`
/// within `tolerance`, or that record no occupancy.
///
/// Every alternate location is kept unless `policy` says otherwise, because the
/// default policy would resolve each group to one conformer before the sums
/// could be taken.
///
/// The value is a dictionary with `intended`, `assessed` and `records`, each
/// record carrying `residue`, `atom_name`, `alternatives`, `assessed`,
/// `occupancy_sum` (`None` when unrecorded) and `issue`
/// (`"missing_occupancy"`, `"sum_mismatch"` or `None`).
#[pyfunction]
#[pyo3(signature = (structure, *, expected_sum, tolerance, policy=None, context=None))]
pub(crate) fn altloc_occupancy_sums(
    py: Python<'_>,
    structure: &PyStructure,
    expected_sum: f64,
    tolerance: f64,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::altloc_occupancy_sums_kernel(AltlocOccupancyOptions {
        expected_sum,
        tolerance,
    });
    let policy = match policy {
        Some(policy) => policy.0.clone(),
        None => molframe::AnalysisPolicy {
            altloc: molframe::AltlocPolicy::KeepAll,
            ..molframe::AnalysisPolicy::default()
        },
    };
    run(
        py,
        structure,
        &policy,
        &kernel,
        |py, report| {
            let value = PyDict::new(py);
            value.set_item("intended", report.intended)?;
            value.set_item("assessed", report.assessed)?;
            let records = PyList::empty(py);
            for record in report.records {
                let entry = PyDict::new(py);
                entry.set_item("residue", record.residue.get())?;
                entry.set_item("atom_name", record.atom_name)?;
                entry.set_item("alternatives", record.alternatives)?;
                entry.set_item("assessed", record.assessed)?;
                entry.set_item("occupancy_sum", record.occupancy_sum)?;
                entry.set_item(
                    "issue",
                    record.issue.map(|issue| match issue {
                        molframe::validation::AltlocOccupancyIssue::MissingOccupancy => {
                            "missing_occupancy"
                        }
                        molframe::validation::AltlocOccupancyIssue::SumMismatch => "sum_mismatch",
                    }),
                )?;
                records.append(entry)?;
            }
            value.set_item("records", records)?;
            Ok(value.into_any().unbind())
        },
        context,
    )
}

/// Bonds of heterogen components that deviate from the covalent expectation
/// by more than `tolerance` ångström.
///
/// The value is a dictionary with `intended`, `assessed` and `outliers`, a
/// table in the layout of `bond_length_deviations`.
#[pyfunction]
#[pyo3(signature = (structure, *, tolerance, policy=None, context=None))]
pub(crate) fn ligand_geometry(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel = molframe::validation::ligand_geometry_kernel(tolerance);
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, report| {
            let value = PyDict::new(py);
            value.set_item("intended", report.intended)?;
            value.set_item("assessed", report.assessed)?;
            let rows = &report.outliers;
            let observed: Vec<f32> = rows.iter().map(|row| row.observed).collect();
            let expected: Vec<f32> = rows.iter().map(|row| row.expected).collect();
            let deviation: Vec<f32> = rows.iter().map(|row| row.deviation).collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("atom_a", &atoms(rows.iter().map(|row| row.atom_a)))
                .indices("atom_b", &atoms(rows.iter().map(|row| row.atom_b)))
                .single("observed", &observed)
                .single("expected", &expected)
                .single("deviation", &deviation)
                .finish();
            value.set_item("outliers", Py::new(py, table)?)?;
            Ok(value.into_any().unbind())
        },
        context,
    )
}

/// The distribution of recorded B factors, and the atoms more than
/// `outlier_standard_deviations` standard deviations from the mean.
///
/// `selection` restricts the atoms assessed. The value is a dictionary with
/// `intended`, `assessed`, `mean`, `variance`, `standard_deviation`,
/// `minimum`, `median`, `maximum` and `outliers` (a table of `atom`, `value`
/// and `z_score`).
#[pyfunction]
#[pyo3(signature = (structure, *, outlier_standard_deviations, selection=None, policy=None, context=None))]
pub(crate) fn b_factor_distribution(
    py: Python<'_>,
    structure: &PyStructure,
    outlier_standard_deviations: f64,
    selection: Option<&crate::bindings::PySelection>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let chosen = match selection {
        Some(selection) => {
            molframe::engine::core::AtomSelection::from_sorted(selection.atom_indices().to_vec())
        }
        None => molframe::engine::core::AtomSelection::All(structure.inner.atom_count()),
    };
    let kernel =
        molframe::validation::b_factor_distribution_kernel(&chosen, outlier_standard_deviations);
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, report| {
            let value = PyDict::new(py);
            value.set_item("intended", report.intended)?;
            value.set_item("assessed", report.assessed)?;
            value.set_item("mean", report.mean)?;
            value.set_item("variance", report.variance)?;
            value.set_item("standard_deviation", report.standard_deviation)?;
            value.set_item("minimum", report.minimum)?;
            value.set_item("median", report.median)?;
            value.set_item("maximum", report.maximum)?;
            let rows = &report.outliers;
            let values: Vec<f64> = rows.iter().map(|row| row.value).collect();
            let scores: Vec<f64> = rows.iter().map(|row| row.z_score).collect();
            let table = TableBuilder::new(py, rows.len())
                .indices("atom", &atoms(rows.iter().map(|row| row.atom)))
                .double("value", &values)
                .double("z_score", &scores)
                .finish();
            value.set_item("outliers", Py::new(py, table)?)?;
            Ok(value.into_any().unbind())
        },
        context,
    )
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(bond_length_deviations, module)?)?;
    module.add_function(wrap_pyfunction!(cis_peptides, module)?)?;
    module.add_function(wrap_pyfunction!(quality_flags, module)?)?;
    module.add_function(wrap_pyfunction!(valence, module)?)?;
    module.add_function(wrap_pyfunction!(planarity, module)?)?;
    module.add_function(wrap_pyfunction!(completeness, module)?)?;
    module.add_function(wrap_pyfunction!(altloc_occupancy_sums, module)?)?;
    module.add_function(wrap_pyfunction!(ligand_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(b_factor_distribution, module)?)
}
