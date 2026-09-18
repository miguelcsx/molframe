//! Native chunk summaries and parent mappings.

use crate::chemistry::PyElement;
use crate::core_topology::PyResidueTable;
use crate::core_values::PyAabb;
use crate::index::PyResidueIndex;
use pyo3::prelude::*;

#[pyclass(name = "ElementMask", from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyElementMask(pub(crate) molframe::core::ElementMask);

#[pyclass(name = "Extremes", from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyExtremes(pub(crate) molframe::core::Extremes);

#[pyclass(name = "AtomChunkStats", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAtomChunkStats(pub(crate) molframe::core::AtomChunkStats);

#[pyclass(name = "ParentMapping", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyParentMapping(pub(crate) molframe::core::ParentMapping);

#[pymethods]
impl PyElementMask {
    #[classattr]
    #[pyo3(name = "EMPTY")]
    fn empty() -> Self {
        Self(molframe::core::ElementMask::EMPTY)
    }

    #[new]
    fn new() -> Self {
        Self::empty()
    }

    fn insert(&mut self, element: PyElement) {
        self.0.insert(element.0);
    }

    fn contains(&self, element: PyElement) -> bool {
        self.0.contains(element.0)
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn union_with(&mut self, other: PyElementMask) {
        self.0.union_with(other.0);
    }

    fn intersects(&self, other: PyElementMask) -> bool {
        self.0.intersects(other.0)
    }

    fn len(&self) -> u32 {
        self.0.len()
    }
}

#[pymethods]
impl PyExtremes {
    #[new]
    fn new() -> Self {
        Self(molframe::core::Extremes::default())
    }

    fn observe(&mut self, value: f32) {
        self.0.observe(value);
    }

    #[getter]
    fn min(&self) -> Option<f32> {
        self.0.min()
    }

    #[getter]
    fn max(&self) -> Option<f32> {
        self.0.max()
    }

    fn excludes_above(&self, at_least: f32) -> bool {
        self.0.excludes_above(at_least)
    }

    fn excludes_below(&self, at_most: f32) -> bool {
        self.0.excludes_below(at_most)
    }
}

#[pymethods]
impl PyAtomChunkStats {
    #[getter]
    fn bounds(&self) -> PyAabb {
        PyAabb(self.0.bounds)
    }

    #[getter]
    fn elements(&self) -> PyElementMask {
        PyElementMask(self.0.elements)
    }

    #[getter]
    fn model(&self) -> u32 {
        self.0.model
    }

    #[getter]
    fn residue_min(&self) -> u32 {
        self.0.residue_min
    }

    #[getter]
    fn residue_max(&self) -> u32 {
        self.0.residue_max
    }

    #[getter]
    fn b_factor(&self) -> PyExtremes {
        PyExtremes(self.0.b_factor)
    }

    #[getter]
    fn occupancy(&self) -> PyExtremes {
        PyExtremes(self.0.occupancy)
    }

    #[getter]
    fn has_altloc(&self) -> bool {
        self.0.has_altloc
    }

    #[getter]
    fn has_missing_coords(&self) -> bool {
        self.0.has_missing_coords
    }

    #[getter]
    fn has_hydrogen(&self) -> bool {
        self.0.has_hydrogen
    }

    fn excludes_element(&self, element: PyElement) -> bool {
        self.0.excludes_element(element.0)
    }

    fn excludes_region(&self, region: PyAabb, distance: f32) -> bool {
        self.0.excludes_region(&region.0, distance)
    }

    fn excludes_residues(&self, first: u32, last: u32) -> bool {
        self.0.excludes_residues(first, last)
    }
}

#[pymethods]
impl PyParentMapping {
    #[staticmethod]
    fn offsets_only() -> Self {
        Self(molframe::core::ParentMapping::OffsetsOnly)
    }

    #[staticmethod]
    fn explicit(parents: Vec<u32>) -> Self {
        Self(molframe::core::ParentMapping::explicit(&parents))
    }

    #[staticmethod]
    fn block_indexed(parents: Vec<u32>) -> Self {
        Self(molframe::core::ParentMapping::block_indexed(&parents))
    }

    #[classattr]
    #[pyo3(name = "DEFAULT_BLOCK")]
    fn default_block() -> u16 {
        molframe::core::ParentMapping::DEFAULT_BLOCK
    }

    fn bytes(&self) -> usize {
        self.0.bytes()
    }

    fn resolve(
        &self,
        local: u32,
        first_atom: u32,
        residues: &PyResidueTable,
    ) -> Option<PyResidueIndex> {
        self.0
            .resolve(local, first_atom, &residues.0)
            .map(PyResidueIndex)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyElementMask>()?;
    module.add_class::<PyExtremes>()?;
    module.add_class::<PyAtomChunkStats>()?;
    module.add_class::<PyParentMapping>()?;
    module.add("TARGET_CHUNK_ATOMS", molframe::core::TARGET_CHUNK_ATOMS)?;
    Ok(())
}
