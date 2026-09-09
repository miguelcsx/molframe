//! Python access to native atom records, chunks, and chunk builders.

use crate::chemistry::PyElement;
use crate::core_chunk_stats::PyAtomChunkStats;
use crate::core_columns::PyPresence;
use crate::core_storage::PyCoordinateBlock;
use crate::core_topology::PyResidueTable;
use crate::core_values::{PyAltId, PyOptionalSymbol, PySymbolId};
use crate::index::PyResidueIndex;
use numpy::ndarray::ArrayView1;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::prelude::*;

#[pyclass(name = "AtomRecord", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAtomRecord {
    position: Option<[f32; 3]>,
    element: PyElement,
    atom_name: PySymbolId,
    auth_atom_name: PyOptionalSymbol,
    alternate_component_id: PyOptionalSymbol,
    alt_id: PyAltId,
    residue: PyResidueIndex,
    occupancy: (f32, PyPresence),
    b_factor: (f32, PyPresence),
    formal_charge: (i8, PyPresence),
    atom_site_id: u32,
}

impl PyAtomRecord {
    fn native(&self) -> pdbiox::AtomRecord {
        pdbiox::AtomRecord {
            position: self.position,
            element: self.element.0,
            atom_name: self.atom_name.0,
            auth_atom_name: self.auth_atom_name.0,
            alternate_component_id: self.alternate_component_id.0,
            alt_id: self.alt_id.0,
            residue: self.residue.0,
            occupancy: (self.occupancy.0, self.occupancy.1.into()),
            b_factor: (self.b_factor.0, self.b_factor.1.into()),
            formal_charge: (self.formal_charge.0, self.formal_charge.1.into()),
            atom_site_id: self.atom_site_id,
        }
    }
}

#[pyclass(name = "AtomChunk", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAtomChunk(pub(crate) pdbiox::AtomChunk);

#[pyclass(name = "ChunkBuilder", skip_from_py_object)]
pub(crate) struct PyChunkBuilder(pub(crate) pdbiox::ChunkBuilder);

fn readonly_f32<'py>(
    py: Python<'py>,
    owner: &PyAtomChunk,
    values: &[f32],
) -> PyResult<Bound<'py, PyArray1<f32>>> {
    let owner = Bound::new(py, owner.clone())?;
    let view = unsafe { ArrayView1::from_shape_ptr(values.len(), values.as_ptr()) };
    // SAFETY: the owner is retained by NumPy, and the native column remains
    // immutable for the lifetime of this cloned chunk.
    let array = unsafe { PyArray1::borrow_from_array(&view, owner.into_any()) };
    array.readwrite().make_nonwriteable();
    Ok(array)
}

#[pymethods]
impl PyAtomRecord {
    #[new]
    #[pyo3(signature = (
        position,
        element,
        atom_name,
        auth_atom_name,
        alternate_component_id,
        alt_id,
        residue,
        occupancy,
        b_factor,
        formal_charge,
        atom_site_id
    ))]
    fn new(
        position: Option<[f32; 3]>,
        element: PyElement,
        atom_name: PySymbolId,
        auth_atom_name: PyOptionalSymbol,
        alternate_component_id: PyOptionalSymbol,
        alt_id: PyAltId,
        residue: PyResidueIndex,
        occupancy: (f32, PyPresence),
        b_factor: (f32, PyPresence),
        formal_charge: (i8, PyPresence),
        atom_site_id: u32,
    ) -> Self {
        Self {
            position,
            element,
            atom_name,
            auth_atom_name,
            alternate_component_id,
            alt_id,
            residue,
            occupancy,
            b_factor,
            formal_charge,
            atom_site_id,
        }
    }

    #[getter]
    fn position(&self) -> Option<[f32; 3]> {
        self.position
    }

    #[getter]
    fn element(&self) -> PyElement {
        self.element.clone()
    }

    #[getter]
    fn atom_name(&self) -> PySymbolId {
        self.atom_name
    }

    #[getter]
    fn auth_atom_name(&self) -> PyOptionalSymbol {
        self.auth_atom_name
    }

    #[getter]
    fn alternate_component_id(&self) -> PyOptionalSymbol {
        self.alternate_component_id
    }

    #[getter]
    fn alt_id(&self) -> PyAltId {
        self.alt_id
    }

    #[getter]
    fn residue(&self) -> PyResidueIndex {
        self.residue
    }

    #[getter]
    fn occupancy(&self) -> (f32, PyPresence) {
        self.occupancy
    }

    #[getter]
    fn b_factor(&self) -> (f32, PyPresence) {
        self.b_factor
    }

    #[getter]
    fn formal_charge(&self) -> (i8, PyPresence) {
        self.formal_charge
    }

    #[getter]
    fn atom_site_id(&self) -> u32 {
        self.atom_site_id
    }
}

#[pymethods]
impl PyAtomChunk {
    fn __len__(&self) -> usize {
        self.0.len() as usize
    }

    fn len(&self) -> usize {
        self.0.len() as usize
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[getter]
    fn atoms(&self) -> [u32; 2] {
        let range = self.0.atoms();
        [range.start, range.end]
    }

    #[getter]
    fn model(&self) -> u32 {
        self.0.model()
    }

    #[getter]
    fn stats(&self) -> PyAtomChunkStats {
        PyAtomChunkStats(*self.0.stats())
    }

    fn positions<'py>(
        &self,
        py: Python<'py>,
        coords: &PyCoordinateBlock,
    ) -> PyResult<Bound<'py, numpy::PyArray2<f32>>> {
        let range = self.0.atoms();
        coords.view(py, Some((range.start, range.end)))
    }

    fn element(&self, local: u32) -> Option<PyElement> {
        self.0.element(local).map(PyElement)
    }

    fn atom_name(&self, local: u32) -> Option<PySymbolId> {
        self.0.atom_name(local).map(PySymbolId)
    }

    fn atom_names_plain(&self) -> Option<Vec<PySymbolId>> {
        self.0
            .atom_names_plain()
            .map(|values| values.iter().copied().map(PySymbolId).collect())
    }

    fn auth_atom_name(&self, local: u32) -> Option<PySymbolId> {
        self.0.auth_atom_name(local).map(PySymbolId)
    }

    fn alternate_component_id(&self, local: u32) -> Option<PySymbolId> {
        self.0.alternate_component_id(local).map(PySymbolId)
    }

    fn alt_id(&self, local: u32) -> Option<PyAltId> {
        self.0.alt_id(local).map(PyAltId)
    }

    fn residue(&self, local: u32, residues: &PyResidueTable) -> Option<PyResidueIndex> {
        self.0.residue(local, &residues.0).map(PyResidueIndex)
    }

    fn occupancy(&self, local: u32) -> Option<(f32, PyPresence)> {
        self.0
            .occupancy(local)
            .map(|(value, presence)| (value, presence.into()))
    }

    fn occupancies_plain<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Option<Bound<'py, PyArray1<f32>>>> {
        self.0
            .occupancies_plain()
            .map(|values| readonly_f32(py, self, values))
            .transpose()
    }

    fn b_factor(&self, local: u32) -> Option<(f32, PyPresence)> {
        self.0
            .b_factor(local)
            .map(|(value, presence)| (value, presence.into()))
    }

    fn b_factors_plain<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f32>>>> {
        self.0
            .b_factors_plain()
            .map(|values| readonly_f32(py, self, values))
            .transpose()
    }

    fn formal_charge(&self, local: u32) -> Option<(i8, PyPresence)> {
        self.0
            .formal_charge(local)
            .map(|(value, presence)| (value, presence.into()))
    }

    fn atom_site_id(&self, local: u32) -> Option<u32> {
        self.0.atom_site_id(local)
    }

    fn has_position(&self, local: u32) -> bool {
        self.0.has_position(local)
    }

    #[pyo3(signature = (local, residues, position=None))]
    fn record(
        &self,
        local: u32,
        residues: &PyResidueTable,
        position: Option<[f32; 3]>,
    ) -> Option<PyAtomRecord> {
        self.0.record(local, &residues.0, position).map(Into::into)
    }
}

#[pymethods]
impl PyChunkBuilder {
    #[new]
    fn new() -> Self {
        Self(pdbiox::ChunkBuilder::new())
    }

    #[staticmethod]
    fn with_target(target: u32) -> Self {
        Self(pdbiox::ChunkBuilder::with_target(target))
    }

    fn reserve(&mut self, atoms: usize) {
        self.0.reserve(atoms);
    }

    fn start_model(&mut self, model: u32) {
        self.0.start_model(model);
    }

    fn push(&mut self, record: PyAtomRecord) {
        self.0.push(record.native());
    }

    fn finish(&mut self, py: Python<'_>) -> (Vec<PyAtomChunk>, PyCoordinateBlock) {
        let builder = std::mem::take(&mut self.0);
        let (chunks, coords) = py.detach(|| builder.finish());
        (
            chunks.into_iter().map(PyAtomChunk).collect(),
            PyCoordinateBlock(coords),
        )
    }
}

impl From<pdbiox::AtomRecord> for PyAtomRecord {
    fn from(value: pdbiox::AtomRecord) -> Self {
        Self {
            position: value.position,
            element: PyElement(value.element),
            atom_name: PySymbolId(value.atom_name),
            auth_atom_name: PyOptionalSymbol(value.auth_atom_name),
            alternate_component_id: PyOptionalSymbol(value.alternate_component_id),
            alt_id: PyAltId(value.alt_id),
            residue: PyResidueIndex(value.residue),
            occupancy: (value.occupancy.0, value.occupancy.1.into()),
            b_factor: (value.b_factor.0, value.b_factor.1.into()),
            formal_charge: (value.formal_charge.0, value.formal_charge.1.into()),
            atom_site_id: value.atom_site_id,
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAtomRecord>()?;
    module.add_class::<PyAtomChunk>()?;
    module.add_class::<PyChunkBuilder>()?;
    Ok(())
}
