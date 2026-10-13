//! Mechanical adapters for the table-backed built-in validation.
//!
//! Every kernel runs in Rust against the bundled reference tables, optionally
//! overlaid by user files that Rust parses. The per-residue results move into
//! numpy arrays without a copy; the report names the tables that were applied.

use crate::bindings::PyStructure;
use crate::execution::{PyExecutionContext, run};
use molframe::validation::{
    BuiltinError, ClashScoreOptions, HydrogenMode, Layout, NucleicPolicy, ReferenceOverrides,
    ReferenceSet, TableInfo, bundled_reference, clash_score, geometry_scores, nucleic_table,
    protein_table,
};
use numpy::{Element, PyArray1, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::{Path, PathBuf};

/// Results of one built-in check: per-residue arrays, counts and provenance.
#[derive(Debug)]
#[pyclass(name = "ValidationReport", frozen, module = "molframe.validation")]
struct PyValidationReport {
    check: &'static str,
    arrays: Py<PyDict>,
    summary: Py<PyDict>,
    tables: Py<PyList>,
}

#[pymethods]
impl PyValidationReport {
    /// Which check produced this report.
    #[getter]
    fn check(&self) -> &'static str {
        self.check
    }

    /// Per-residue (or per-term) numpy arrays by name, shared with the report.
    #[getter]
    fn arrays<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        self.arrays.bind(py).clone()
    }

    /// Counts and rates folded from the arrays.
    #[getter]
    fn summary<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        self.summary.bind(py).clone()
    }

    /// The tables or policy applied: name, version, source and licence.
    #[getter]
    fn tables<'py>(&self, py: Python<'py>) -> Bound<'py, PyList> {
        self.tables.bind(py).clone()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        match self.arrays.bind(py).get_item(name)? {
            Some(array) => Ok(array),
            None => Err(crate::error::key(format!("no array named `{name}`"))),
        }
    }

    fn __repr__(&self) -> String {
        format!("ValidationReport(check={:?})", self.check)
    }
}

fn table_error(error: &BuiltinError) -> PyErr {
    crate::error::value(error.to_string())
}

fn array<T: Element>(py: Python<'_>, values: Vec<T>) -> Bound<'_, PyArray1<T>> {
    PyArray1::from_vec(py, values)
}

fn matrix(py: Python<'_>, values: Vec<f32>, columns: usize) -> PyResult<Bound<'_, PyAny>> {
    let rows = values
        .len()
        .checked_div(columns)
        .ok_or_else(|| crate::error::value("matrix needs columns"))?;
    Ok(array(py, values).reshape([rows, columns])?.into_any())
}

fn info_list<'py>(py: Python<'py>, tables: &[TableInfo]) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty(py);
    for info in tables {
        let entry = PyDict::new(py);
        entry.set_item("name", &info.name)?;
        entry.set_item("version", &info.version)?;
        entry.set_item("source", &info.source)?;
        entry.set_item("licence", &info.licence)?;
        list.append(entry)?;
    }
    Ok(list)
}

fn counts<'py, V: IntoPyObject<'py>>(
    py: Python<'py>,
    pairs: Vec<(&str, V)>,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for (name, value) in pairs {
        dict.set_item(name, value)?;
    }
    Ok(dict)
}

fn report(
    check: &'static str,
    arrays: Bound<'_, PyDict>,
    summary: Bound<'_, PyDict>,
    tables: Bound<'_, PyList>,
) -> PyValidationReport {
    PyValidationReport {
        check,
        arrays: arrays.unbind(),
        summary: summary.unbind(),
        tables: tables.unbind(),
    }
}

fn reference(
    rama_table: Option<&Path>,
    rotamer_table: Option<&Path>,
    geometry_table: Option<&Path>,
    chi_table: Option<&Path>,
) -> PyResult<ReferenceSet> {
    let overrides =
        ReferenceOverrides::from_paths(rama_table, rotamer_table, geometry_table, chi_table)
            .map_err(|error| table_error(&error))?;
    let base = bundled_reference().map_err(|error| table_error(&error))?;
    Ok(base.clone().overlaid(overrides))
}

/// Ramachandran, rotamer and omega classes for every residue.
///
/// Reference tables default to the bundled Top8000 grids; each `*_table`
/// path replaces one table with a user file parsed in Rust.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
#[pyo3(signature = (structure, *, rama_table=None, rotamer_table=None, chi_table=None, context=None))]
fn protein_conformation(
    py: Python<'_>,
    structure: &PyStructure,
    rama_table: Option<PathBuf>,
    rotamer_table: Option<PathBuf>,
    chi_table: Option<PathBuf>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyValidationReport> {
    let set = reference(
        rama_table.as_deref(),
        rotamer_table.as_deref(),
        None,
        chi_table.as_deref(),
    )?;
    let engine = structure.inner.engine();
    let table = run(py, context, |ctx| {
        protein_table(&Layout::new(engine), &set, ctx)
    })?
    .map_err(|error| table_error(&error))?;
    let s = table.summary();
    let arrays = PyDict::new(py);
    arrays.set_item("phi", array(py, table.phi))?;
    arrays.set_item("psi", array(py, table.psi))?;
    arrays.set_item("omega", array(py, table.omega))?;
    arrays.set_item("chi", matrix(py, table.chi, 4)?)?;
    arrays.set_item("rama_category", array(py, table.rama_category))?;
    arrays.set_item("rama_density", array(py, table.rama_density))?;
    arrays.set_item("rama_class", array(py, table.rama_class))?;
    arrays.set_item("rotamer_class", array(py, table.rotamer_class))?;
    arrays.set_item("omega_class", array(py, table.omega_class))?;
    let summary = counts(
        py,
        vec![
            ("rama_evaluated", s.rama_evaluated),
            ("rama_favoured", s.rama_favoured),
            ("rama_allowed", s.rama_allowed),
            ("rama_outliers", s.rama_outliers),
            ("rotamer_evaluated", s.rotamer_evaluated),
            ("rotamer_favoured", s.rotamer_favoured),
            ("rotamer_allowed", s.rotamer_allowed),
            ("rotamer_outliers", s.rotamer_outliers),
            ("omega_evaluated", s.omega_evaluated),
            ("cis_peptides", s.cis_peptides),
            ("twisted_peptides", s.twisted_peptides),
        ],
    )?;
    let tables = info_list(py, &table.tables)?;
    Ok(report("protein", arrays, summary, tables))
}

/// Bond-length and bond-angle z-scores against Engh and Huber targets.
///
/// `cutoff` selects the terms listed in the `outlier_*` arrays.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
#[pyo3(signature = (structure, *, cutoff=4.0, geometry_table=None, chi_table=None, context=None))]
fn geometry_z(
    py: Python<'_>,
    structure: &PyStructure,
    cutoff: f64,
    geometry_table: Option<PathBuf>,
    chi_table: Option<PathBuf>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyValidationReport> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(crate::error::value("cutoff must be finite and positive"));
    }
    let set = reference(None, None, geometry_table.as_deref(), chi_table.as_deref())?;
    let engine = structure.inner.engine();
    let scores = run(py, context, |ctx| {
        geometry_scores(&Layout::new(engine), &set, cutoff, ctx)
    })?
    .map_err(|error| table_error(&error))?;
    let terms = &scores.outliers;
    let arrays = PyDict::new(py);
    arrays.set_item(
        "outlier_residue",
        array(py, terms.iter().map(|t| t.residue).collect()),
    )?;
    arrays.set_item("outlier_z", array(py, terms.iter().map(|t| t.z).collect()))?;
    arrays.set_item(
        "outlier_observed",
        array(py, terms.iter().map(|t| t.observed).collect()),
    )?;
    arrays.set_item(
        "outlier_expected",
        array(py, terms.iter().map(|t| t.expected).collect()),
    )?;
    arrays.set_item(
        "outlier_is_angle",
        array(py, terms.iter().map(|t| u8::from(t.is_angle())).collect()),
    )?;
    let summary = PyDict::new(py);
    summary.set_item("bond_count", scores.bond_count)?;
    summary.set_item("angle_count", scores.angle_count)?;
    summary.set_item("bond_rmsz", scores.bond_rmsz)?;
    summary.set_item("angle_rmsz", scores.angle_rmsz)?;
    summary.set_item("outliers", terms.len())?;
    summary.set_item("cutoff", cutoff)?;
    let tables = info_list(py, &scores.tables)?;
    arrays.set_item("worst_bond_z", array(py, scores.worst_bond_z))?;
    arrays.set_item("worst_angle_z", array(py, scores.worst_angle_z))?;
    Ok(report("geometry-z", arrays, summary, tables))
}

/// Sugar pucker and backbone sanity per nucleotide.
///
/// Suite sanity flags delta/pucker mismatch and gamma outside the staggered
/// windows; it is a filter, not a named suite clustering.
#[pyfunction]
#[pyo3(signature = (structure, *, context=None))]
fn nucleic_conformation(
    py: Python<'_>,
    structure: &PyStructure,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyValidationReport> {
    let engine = structure.inner.engine();
    let table = run(py, context, |ctx| {
        nucleic_table(&Layout::new(engine), &NucleicPolicy::default(), ctx)
    })?
    .map_err(|error| table_error(&error))?;
    let s = table.summary();
    let arrays = PyDict::new(py);
    arrays.set_item("nu", matrix(py, table.nu, 5)?)?;
    arrays.set_item("pucker_phase", array(py, table.pucker_phase))?;
    arrays.set_item("pucker_amplitude", array(py, table.pucker_amplitude))?;
    arrays.set_item("pucker_class", array(py, table.pucker_class))?;
    let columns = molframe::validation::NUCLEIC_TORSIONS;
    arrays.set_item("torsions", matrix(py, table.torsions, columns)?)?;
    arrays.set_item("suite_flags", array(py, table.suite_flags))?;
    arrays.set_item("suite_class", array(py, table.suite_class))?;
    let summary = counts(
        py,
        vec![
            ("pucker_evaluated", s.pucker_evaluated),
            ("c3_endo", s.c3_endo),
            ("c2_endo", s.c2_endo),
            ("other_pucker", s.other_pucker),
            ("suite_evaluated", s.suite_evaluated),
            ("suite_flagged", s.suite_flagged),
        ],
    )?;
    let tables = info_list(py, &table.tables)?;
    Ok(report("nucleic-conformation", arrays, summary, tables))
}

/// Clashes per 1000 atoms, heavy-atom or hydrogen-aware.
///
/// `hydrogens="auto"` uses explicit hydrogens when the model has any;
/// `"ignore"` forces the heavy-atom reading. The summary records the policy.
#[pyfunction]
#[pyo3(signature = (structure, *, tolerance=0.4, radii="bondi", hydrogens="auto", polar_margin=0.6, context=None))]
fn clashscore(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    radii: &str,
    hydrogens: &str,
    polar_margin: f32,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyValidationReport> {
    let mode = match hydrogens {
        "auto" => HydrogenMode::Auto,
        "ignore" => HydrogenMode::Ignore,
        _ => return Err(crate::error::value("hydrogens must be `auto` or `ignore`")),
    };
    let options = ClashScoreOptions {
        tolerance,
        radius_set: crate::chemistry::radius_set(radii)?,
        hydrogens: mode,
        polar_margin,
    };
    let engine = structure.inner.engine();
    let score = run(py, context, |ctx| clash_score(engine, options, ctx))?
        .map_err(|error| table_error(&error))?;
    let arrays = PyDict::new(py);
    arrays.set_item("residue_clashes", array(py, score.residue_clashes))?;
    let summary = PyDict::new(py);
    summary.set_item("score", score.score)?;
    summary.set_item("clashes", score.clashes)?;
    summary.set_item("atoms", score.atoms)?;
    summary.set_item("worst_overlap", score.worst_overlap)?;
    summary.set_item("hydrogen_aware", score.hydrogen_aware)?;
    summary.set_item("tolerance", score.options.tolerance)?;
    summary.set_item("polar_margin", score.options.polar_margin)?;
    let policy = TableInfo::new(
        "clashscore-policy",
        "molframe-1",
        format!(
            "radii={radii} tolerance={} polar_margin={} hydrogen_aware={}",
            score.options.tolerance, score.options.polar_margin, score.hydrogen_aware
        ),
        "not applicable",
    );
    let tables = info_list(py, &[policy])?;
    Ok(report("clashscore", arrays, summary, tables))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyValidationReport>()?;
    module.add_function(wrap_pyfunction!(protein_conformation, module)?)?;
    module.add_function(wrap_pyfunction!(geometry_z, module)?)?;
    module.add_function(wrap_pyfunction!(nucleic_conformation, module)?)?;
    module.add_function(wrap_pyfunction!(clashscore, module)?)
}
