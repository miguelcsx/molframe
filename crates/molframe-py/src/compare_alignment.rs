//! Structure alignment by combinatorial extension, contact-area difference scores, contact
//! map overlap and the sequence-based mapping of chains and residues.

use crate::execution::PyExecutionContext;
use crate::table::{PyTable, TableBuilder};
use molframe::compare::{
    CeOptions, CeSignificanceProfile, ContactArea, ce_alignments, contact_map_similarity,
};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1};
use pyo3::prelude::*;

fn coordinates<'a>(array: &'a numpy::PyReadonlyArray2<'_, f32>) -> PyResult<&'a [[f32; 3]]> {
    crate::bindings::coordinates(array)
}

/// One combinatorial-extension correspondence between two sets of guide atoms.
#[derive(Clone, Debug)]
#[pyclass(
    name = "CeAlignment",
    frozen,
    skip_from_py_object,
    module = "molframe.compare"
)]
pub(crate) struct PyCeAlignment {
    inner: molframe::compare::CeAlignment,
}

#[pymethods]
impl PyCeAlignment {
    /// Indices into the reference guide atoms.
    #[getter]
    fn reference_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        let indices: Vec<u64> = self
            .inner
            .reference_indices
            .iter()
            .map(|&i| i as u64)
            .collect();
        indices.into_pyarray(py)
    }

    /// The matching indices into the mobile guide atoms.
    #[getter]
    fn mobile_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        let indices: Vec<u64> = self
            .inner
            .mobile_indices
            .iter()
            .map(|&i| i as u64)
            .collect();
        indices.into_pyarray(py)
    }

    /// Aligned fragment pairs before they are expanded to guide atoms.
    #[getter]
    const fn fragment_count(&self) -> usize {
        self.inner.fragment_count
    }

    /// The path similarity; values nearer zero are better.
    #[getter]
    const fn similarity(&self) -> f64 {
        self.inner.similarity
    }

    /// The empirical significance estimate, when a calibration profile was chosen.
    #[getter]
    const fn z_score(&self) -> Option<f64> {
        self.inner.z_score
    }

    /// RMSD of the correspondence after one shared rigid fit.
    #[getter]
    const fn rmsd(&self) -> f64 {
        self.inner.rmsd
    }

    fn __len__(&self) -> usize {
        self.inner.reference_indices.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "CeAlignment(aligned={}, rmsd={:.3})",
            self.inner.reference_indices.len(),
            self.inner.rmsd
        )
    }
}

/// Combinatorial-extension alignments of two sets of guide atoms (for instance C-alpha),
/// best first: the most atoms aligned, then the least RMSD.
///
/// The defaults are the original CE search: eight-atom fragments, gaps up to 30, 20 paths,
/// fragment and path similarity thresholds of -3 and -4, and the original window-eight
/// significance calibration. `significance=False` reports no z-score. `memory_limit` bounds the
/// search workspace in bytes.
#[pyfunction]
#[pyo3(signature = (
    reference,
    mobile,
    *,
    window_size=8,
    max_gap=30,
    max_paths=20,
    fragment_threshold=-3.0,
    path_threshold=-4.0,
    significance=true,
    memory_limit=100_000_000,
))]
#[allow(clippy::too_many_arguments)]
fn ce_align(
    py: Python<'_>,
    reference: &Bound<'_, PyArray2<f32>>,
    mobile: &Bound<'_, PyArray2<f32>>,
    window_size: usize,
    max_gap: usize,
    max_paths: usize,
    fragment_threshold: f64,
    path_threshold: f64,
    significance: bool,
    memory_limit: usize,
) -> PyResult<Vec<PyCeAlignment>> {
    let (reference, mobile) = (reference.readonly(), mobile.readonly());
    let (reference, mobile) = (coordinates(&reference)?, coordinates(&mobile)?);
    let options = CeOptions {
        window_size,
        max_gap,
        max_paths,
        fragment_similarity_threshold: fragment_threshold,
        path_similarity_threshold: path_threshold,
        significance: significance.then_some(CeSignificanceProfile::OriginalWindowEight),
        memory_limit_bytes: memory_limit,
    };
    let found = py
        .detach(|| ce_alignments(reference, mobile, options))
        .map_err(crate::error::kernel)?;
    Ok(found
        .into_iter()
        .map(|inner| PyCeAlignment { inner })
        .collect())
}

fn column<'py, T: numpy::Element + Copy>(
    table: &Bound<'py, PyAny>,
    name: &str,
) -> PyResult<Vec<T>> {
    let array: PyReadonlyArray1<'py, T> = table.get_item(name)?.extract()?;
    Ok(array.as_array().to_vec())
}

fn areas(table: &Bound<'_, PyAny>) -> PyResult<Vec<ContactArea>> {
    let first: Vec<u32> = column(table, "first")?;
    let second: Vec<u32> = column(table, "second")?;
    let area: Vec<f64> = column(table, "area")?;
    if first.len() != second.len() || first.len() != area.len() {
        return Err(crate::error::value(
            "first, second and area must have the same length",
        ));
    }
    Ok(first
        .iter()
        .zip(&second)
        .zip(&area)
        .map(|((&first, &second), &area)| ContactArea {
            first,
            second,
            area,
        })
        .collect())
}

/// Contact areas between residues: the surface area on each atom's solvent-expanded
/// sphere that lies inside a neighbouring atom's, summed by residue pair.
///
/// `residues` gives each atom's residue identifier (an integer you choose so that two
/// structures to be compared use the same ones). Columns: `first`, `second`, `area`.
#[pyfunction]
#[pyo3(signature = (positions, radii, residues, *, probe=1.4, density=4.0, context=None))]
fn contact_areas(
    py: Python<'_>,
    positions: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    residues: &Bound<'_, PyArray1<u32>>,
    probe: f32,
    density: f32,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyTable> {
    let positions = positions.readonly();
    let radii = radii.readonly();
    let residues = residues.readonly();
    let positions = coordinates(&positions)?;
    let (radii, residues) = (
        radii
            .as_slice()
            .map_err(|_| crate::error::value("radii must be C-contiguous"))?,
        residues
            .as_slice()
            .map_err(|_| crate::error::value("residues must be C-contiguous"))?,
    );
    let found = crate::execution::run(py, context, |context| {
        molframe::compare::cad_contact_areas(positions, radii, residues, probe, density, context)
    })?
    .map_err(crate::error::kernel)?;
    let first: Vec<u32> = found.iter().map(|contact| contact.first).collect();
    let second: Vec<u32> = found.iter().map(|contact| contact.second).collect();
    let area: Vec<f64> = found.iter().map(|contact| contact.area).collect();
    Ok(TableBuilder::new(py, found.len())
        .indices("first", &first)
        .indices("second", &second)
        .double("area", &area)
        .finish())
}

/// The contact-area difference score of a model against a reference.
#[derive(Clone, Debug)]
#[pyclass(
    name = "CadScore",
    frozen,
    skip_from_py_object,
    module = "molframe.compare"
)]
pub(crate) struct PyCadScore {
    inner: molframe::compare::CadScore,
}

#[pymethods]
impl PyCadScore {
    /// The global score in `[0, 1]`.
    #[getter]
    const fn score(&self) -> f64 {
        self.inner.score
    }

    /// The reference area, the global denominator.
    #[getter]
    const fn reference_area(&self) -> f64 {
        self.inner.reference_area
    }

    /// The capped area difference, the global numerator.
    #[getter]
    const fn lost_area(&self) -> f64 {
        self.inner.lost_area
    }

    /// Every reference contact and its contribution. Columns: `first`, `second`,
    /// `reference_area`, `model_area`, `lost_area`.
    fn contacts(&self, py: Python<'_>) -> PyTable {
        let rows = &self.inner.contacts;
        let first: Vec<u32> = rows.iter().map(|row| row.first).collect();
        let second: Vec<u32> = rows.iter().map(|row| row.second).collect();
        let reference: Vec<f64> = rows.iter().map(|row| row.reference_area).collect();
        let model: Vec<f64> = rows.iter().map(|row| row.model_area).collect();
        let lost: Vec<f64> = rows.iter().map(|row| row.lost_area).collect();
        TableBuilder::new(py, rows.len())
            .indices("first", &first)
            .indices("second", &second)
            .double("reference_area", &reference)
            .double("model_area", &model)
            .double("lost_area", &lost)
            .finish()
    }

    /// The score of each residue. Columns: `residue`, `reference_area`, `lost_area`, `score`.
    fn local(&self, py: Python<'_>) -> PyTable {
        let rows = &self.inner.local;
        let residue: Vec<u32> = rows.iter().map(|row| row.residue).collect();
        let reference: Vec<f64> = rows.iter().map(|row| row.reference_area).collect();
        let lost: Vec<f64> = rows.iter().map(|row| row.lost_area).collect();
        let score: Vec<f64> = rows.iter().map(|row| row.score).collect();
        TableBuilder::new(py, rows.len())
            .indices("residue", &residue)
            .double("reference_area", &reference)
            .double("lost_area", &lost)
            .double("score", &score)
            .finish()
    }
}

/// CAD score of a model's contact areas against a reference's.
///
/// Each argument has columns `first`, `second` and `area` (as `contact_areas` returns),
/// with residue identifiers that mean the same residue in both.
#[pyfunction]
fn cad_score(
    py: Python<'_>,
    reference: &Bound<'_, PyAny>,
    model: &Bound<'_, PyAny>,
) -> PyResult<PyCadScore> {
    let (reference, model) = (areas(reference)?, areas(model)?);
    py.detach(|| molframe::compare::cad_score(&reference, &model))
        .map(|inner| PyCadScore { inner })
        .map_err(crate::error::kernel)
}

/// Overlap of two residue contact maps.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "ContactSimilarity",
    frozen,
    skip_from_py_object,
    module = "molframe.compare"
)]
pub(crate) struct PyContactSimilarity {
    inner: molframe::compare::ContactSimilarity,
}

#[pymethods]
impl PyContactSimilarity {
    /// Contacts present in both maps.
    #[getter]
    const fn shared(&self) -> usize {
        self.inner.shared
    }

    /// Contacts present in either.
    #[getter]
    const fn union(&self) -> usize {
        self.inner.union
    }

    /// `shared / union`; one when both maps are empty.
    #[getter]
    const fn jaccard(&self) -> f64 {
        self.inner.jaccard
    }
}

/// Jaccard overlap of two contact maps, each an `(n, 2)` array of residue pairs; the order
/// within a pair does not matter.
#[pyfunction]
fn contact_similarity(
    py: Python<'_>,
    first: &Bound<'_, PyArray2<u32>>,
    second: &Bound<'_, PyArray2<u32>>,
) -> PyResult<PyContactSimilarity> {
    let pairs = |array: &Bound<'_, PyArray2<u32>>| -> PyResult<Vec<(u32, u32)>> {
        let readonly = array.readonly();
        let view = readonly.as_array();
        if view.ncols() != 2 {
            return Err(crate::error::value("a contact map must have two columns"));
        }
        Ok(view
            .rows()
            .into_iter()
            .map(|row| (row[0], row[1]))
            .collect())
    };
    let (first, second) = (pairs(first)?, pairs(second)?);
    Ok(PyContactSimilarity {
        inner: py.detach(|| contact_map_similarity(&first, &second)),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCeAlignment>()?;
    module.add_class::<PyCadScore>()?;
    module.add_class::<PyContactSimilarity>()?;
    module.add_function(wrap_pyfunction!(ce_align, module)?)?;
    module.add_function(wrap_pyfunction!(contact_areas, module)?)?;
    module.add_function(wrap_pyfunction!(cad_score, module)?)?;
    module.add_function(wrap_pyfunction!(contact_similarity, module)?)?;
    crate::compare_mapping::register(module)
}
