//! Owned semantic structure-difference projections.

use crate::core_contract::PyStructureDifferenceOptions;
use crate::structure::PyStructure;
use pyo3::prelude::*;

type CellValues = ([f64; 3], [f64; 3]);
type OptionalCellChange = (Option<CellValues>, Option<CellValues>);

#[pyclass(name = "CountDifference", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCountDifference {
    #[pyo3(get)]
    pub(crate) left: usize,
    #[pyo3(get)]
    pub(crate) right: usize,
}

#[pyclass(name = "MetadataDifference", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMetadataDifference {
    #[pyo3(get)]
    pub(crate) id: Option<(Option<String>, Option<String>)>,
    #[pyo3(get)]
    pub(crate) title: Option<(Option<String>, Option<String>)>,
    #[pyo3(get)]
    pub(crate) method: Option<(Option<String>, Option<String>)>,
    #[pyo3(get)]
    pub(crate) resolution: Option<(Option<f32>, Option<f32>)>,
    #[pyo3(get)]
    pub(crate) cell: Option<OptionalCellChange>,
}

#[pyclass(name = "StructureDifference", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructureDifference {
    #[pyo3(get)]
    pub(crate) metadata: PyMetadataDifference,
    #[pyo3(get)]
    pub(crate) models: Option<PyCountDifference>,
    #[pyo3(get)]
    pub(crate) entities: Option<PyCountDifference>,
    #[pyo3(get)]
    pub(crate) chains: Option<PyCountDifference>,
    #[pyo3(get)]
    pub(crate) residues: Option<PyCountDifference>,
    #[pyo3(get)]
    pub(crate) atoms: Option<PyCountDifference>,
    #[pyo3(get)]
    pub(crate) bonds: Option<PyCountDifference>,
    #[pyo3(get)]
    pub(crate) changed_entities: usize,
    #[pyo3(get)]
    pub(crate) changed_chains: usize,
    #[pyo3(get)]
    pub(crate) changed_residues: usize,
    #[pyo3(get)]
    pub(crate) changed_atoms: usize,
    #[pyo3(get)]
    pub(crate) changed_bonds: usize,
    #[pyo3(get)]
    pub(crate) changed_positions: usize,
    #[pyo3(get)]
    pub(crate) maximum_displacement: Option<f32>,
}

#[pymethods]
impl PyStructureDifference {
    fn is_empty(&self) -> bool {
        self.metadata.id.is_none()
            && self.metadata.title.is_none()
            && self.metadata.method.is_none()
            && self.metadata.resolution.is_none()
            && self.metadata.cell.is_none()
            && self.models.is_none()
            && self.entities.is_none()
            && self.chains.is_none()
            && self.residues.is_none()
            && self.atoms.is_none()
            && self.bonds.is_none()
            && self.changed_entities == 0
            && self.changed_chains == 0
            && self.changed_residues == 0
            && self.changed_atoms == 0
            && self.changed_bonds == 0
            && self.changed_positions == 0
    }
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (other, coordinate_tolerance=0.0, options=None))]
    fn difference(
        &self,
        other: &PyStructure,
        coordinate_tolerance: f32,
        options: Option<&PyStructureDifferenceOptions>,
    ) -> PyResult<PyStructureDifference> {
        structure_difference_with_options(self, other, coordinate_tolerance, options)
    }
}

#[pyfunction]
#[pyo3(signature = (left, right, coordinate_tolerance=0.0, options=None))]
pub(crate) fn structure_difference(
    py: Python<'_>,
    left: &PyStructure,
    right: &PyStructure,
    coordinate_tolerance: f32,
    options: Option<&PyStructureDifferenceOptions>,
) -> PyResult<PyStructureDifference> {
    py.detach(move || -> PyResult<PyStructureDifference> {
        structure_difference_with_options(left, right, coordinate_tolerance, options)
    })
}

fn structure_difference_with_options(
    left: &PyStructure,
    right: &PyStructure,
    coordinate_tolerance: f32,
    options: Option<&PyStructureDifferenceOptions>,
) -> PyResult<PyStructureDifference> {
    let coordinate_tolerance =
        options.map_or(coordinate_tolerance, |value| value.0.coordinate_tolerance);
    pdbiox::structure_difference(
        left.structure(),
        right.structure(),
        pdbiox::StructureDifferenceOptions {
            coordinate_tolerance,
        },
    )
    .map(Into::into)
    .map_err(|error| crate::errors::DifferenceError::new_err(error.to_string()))
}

impl From<pdbiox::StructureDifference> for PyStructureDifference {
    fn from(value: pdbiox::StructureDifference) -> Self {
        Self {
            metadata: value.metadata.into(),
            models: value.models.map(Into::into),
            entities: value.entities.map(Into::into),
            chains: value.chains.map(Into::into),
            residues: value.residues.map(Into::into),
            atoms: value.atoms.map(Into::into),
            bonds: value.bonds.map(Into::into),
            changed_entities: value.changed_entities,
            changed_chains: value.changed_chains,
            changed_residues: value.changed_residues,
            changed_atoms: value.changed_atoms,
            changed_bonds: value.changed_bonds,
            changed_positions: value.changed_positions,
            maximum_displacement: value.maximum_displacement,
        }
    }
}

impl From<pdbiox::CountDifference> for PyCountDifference {
    fn from(value: pdbiox::CountDifference) -> Self {
        Self {
            left: value.left,
            right: value.right,
        }
    }
}

impl From<pdbiox::MetadataDifference> for PyMetadataDifference {
    fn from(value: pdbiox::MetadataDifference) -> Self {
        Self {
            id: text_change(value.id),
            title: text_change(value.title),
            method: text_change(value.method),
            resolution: value.resolution.map(|change| (change.left, change.right)),
            cell: value.cell.map(|change| {
                (
                    change.left.map(|cell| (cell.lengths, cell.angles)),
                    change.right.map(|cell| (cell.lengths, cell.angles)),
                )
            }),
        }
    }
}

fn text_change(
    value: Option<pdbiox::ValueDifference<Option<Box<str>>>>,
) -> Option<(Option<String>, Option<String>)> {
    value.map(|change| {
        (
            change.left.map(|value| value.to_string()),
            change.right.map(|value| value.to_string()),
        )
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCountDifference>()?;
    module.add_class::<PyMetadataDifference>()?;
    module.add_class::<PyStructureDifference>()?;
    module.add_function(wrap_pyfunction!(structure_difference, module)?)?;
    Ok(())
}
