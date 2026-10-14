//! Mechanical adapters for the table-backed built-in validation.
//!
//! Every kernel runs in Rust against reference tables that Rust resolves: a
//! `*_table` path, else the profile's entry in the cache, else a coded error that
//! says how to supply one. No table is built in. The per-residue results move into
//! numpy arrays without a copy; the report names the tables that were applied.

use crate::bindings::PyStructure;
use crate::execution::{PyExecutionContext, run};
use molframe::chemistry::{TableIdentity, TableRequest};
use molframe::validation::{
    ClashScoreOptions, HydrogenMode, Layout, ReferenceSet, clash_score, geometry_scores,
    nucleic_table as score_nucleic, protein_table,
};
use numpy::{Element, PyArray1, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::PathBuf;

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

fn info_list<'py>(py: Python<'py>, tables: &[TableIdentity]) -> PyResult<Bound<'py, PyList>> {
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

fn named<T>(path: Option<&PathBuf>) -> TableRequest<T> {
    TableRequest::choose(path.map(|path| path.to_string_lossy()).as_deref())
}

struct Paths {
    rama: Option<PathBuf>,
    rotamer: Option<PathBuf>,
    geometry: Option<PathBuf>,
    chi: Option<PathBuf>,
    nucleic: Option<PathBuf>,
}

fn reference(paths: &Paths) -> ReferenceSet {
    let mut set = ReferenceSet::default();
    set.rama = named(paths.rama.as_ref());
    set.rotamers = named(paths.rotamer.as_ref());
    set.geometry = named(paths.geometry.as_ref());
    set.chi = named(paths.chi.as_ref());
    set.nucleic = named(paths.nucleic.as_ref());
    set
}

/// Ramachandran, rotamer and omega classes for every residue.
///
/// Each `*_table` path names a table file parsed in Rust; without one the profile's
/// entry in the table cache is used, and a missing table raises a coded error.
#[pyfunction]
#[pyo3(signature = (structure, *, rama_table=None, rotamer_table=None, chi_table=None, context=None))]
fn protein_conformation(
    py: Python<'_>,
    structure: &PyStructure,
    rama_table: Option<PathBuf>,
    rotamer_table: Option<PathBuf>,
    chi_table: Option<PathBuf>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyValidationReport> {
    let set = reference(&Paths {
        rama: rama_table,
        rotamer: rotamer_table,
        geometry: None,
        chi: chi_table,
        nucleic: None,
    });
    let engine = structure.inner.engine();
    let table = run(py, context, |ctx| {
        protein_table(&Layout::new(engine), &set, ctx)
    })?
    .map_err(crate::error::kernel)?;
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

/// Bond-length and bond-angle z-scores against the selected geometry table.
///
/// `cutoff` selects the terms listed in the `outlier_*` arrays.
#[pyfunction]
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
    let set = reference(&Paths {
        rama: None,
        rotamer: None,
        geometry: geometry_table,
        chi: chi_table,
        nucleic: None,
    });
    let engine = structure.inner.engine();
    let scores = run(py, context, |ctx| {
        geometry_scores(&Layout::new(engine), &set, cutoff, ctx)
    })?
    .map_err(crate::error::kernel)?;
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
#[pyo3(signature = (structure, *, nucleic_table=None, context=None))]
fn nucleic_conformation(
    py: Python<'_>,
    structure: &PyStructure,
    nucleic_table: Option<PathBuf>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyValidationReport> {
    let set = reference(&Paths {
        rama: None,
        rotamer: None,
        geometry: None,
        chi: None,
        nucleic: nucleic_table,
    });
    let engine = structure.inner.engine();
    let table = run(py, context, |ctx| {
        score_nucleic(&Layout::new(engine), &set, ctx)
    })?
    .map_err(crate::error::kernel)?;
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
#[pyo3(signature = (structure, *, tolerance=0.4, radii=None, hydrogens="auto", polar_margin=0.6, context=None))]
fn clashscore(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    radii: Option<&str>,
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
        radii: molframe::chemistry::TableRequest::<molframe::chemistry::VdwRadii>::choose(radii),
        hydrogens: mode,
        polar_margin,
    };
    let engine = structure.inner.engine();
    let store = molframe::chemistry::TableStore::global();
    let result = run(py, context, |ctx| {
        clash_score(engine, &options, &store, ctx)
    })?
    .map_err(crate::error::kernel)?;
    let arrays = PyDict::new(py);
    arrays.set_item("residue_clashes", array(py, result.residue_clashes))?;
    let summary = PyDict::new(py);
    summary.set_item("score", result.score)?;
    summary.set_item("clashes", result.clashes)?;
    summary.set_item("atoms", result.atoms)?;
    summary.set_item("worst_overlap", result.worst_overlap)?;
    summary.set_item("hydrogen_aware", result.hydrogen_aware)?;
    summary.set_item("tolerance", result.options.tolerance)?;
    summary.set_item("polar_margin", result.options.polar_margin)?;
    let policy = TableIdentity::supplied(
        "clashscore-policy",
        "molframe-1",
        &format!(
            "radii={} tolerance={} polar_margin={} hydrogen_aware={}",
            result.radii,
            result.options.tolerance,
            result.options.polar_margin,
            result.hydrogen_aware
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
