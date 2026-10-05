//! Small molecules as MDL records: read, build, write and convert to structures.

use crate::bindings::PyStructure;
use molframe::Element;
use molframe::chemistry::{
    MolAtom, MolAtomMetadata, MolBond, MolBondMetadata, MolRecord, MolVersion, Molecule,
    SdfProperty, SmartsPattern, mol_record_to_structure, parse_mol_record, parse_sdf_records,
    structure_to_molecule, write_mol, write_sdf,
};
use molframe::{Code, Diagnostic};
use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedBytes;
use pyo3::types::{PyDict, PyDictMethods};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

/// The text of a path, a string of record text, or bytes.
///
/// A `str` that names an existing file is read as a path; any other `str` is the
/// record text itself, so a literal MOL block and a file name are both accepted.
fn text_of(source: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(bytes) = source.extract::<PyBackedBytes>() {
        return String::from_utf8(bytes.to_vec())
            .map_err(|_| crate::error::value("the record is not UTF-8 text"));
    }
    if let Ok(path) = source.extract::<PathBuf>()
        && path.is_file()
    {
        return std::fs::read_to_string(&path).map_err(|error| {
            crate::error::kernel(
                Diagnostic::new(Code::E7101).with_message(format!("{}: {error}", path.display())),
            )
        });
    }
    match source.extract::<String>() {
        Ok(text) => Ok(text),
        Err(_) => Err(crate::error::type_error(
            "source must be a path, text or bytes",
        )),
    }
}

/// One MOL block or SDF record: a molecular graph, its header and its data fields.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Molecule",
    frozen,
    skip_from_py_object,
    module = "molframe.chemistry"
)]
pub(crate) struct PyMolecule {
    record: Arc<MolRecord>,
}

impl PyMolecule {
    fn new(record: MolRecord) -> Self {
        Self {
            record: Arc::new(record),
        }
    }
}

fn element_of(symbol: &str) -> PyResult<Element> {
    Element::from_symbol(symbol)
        .ok_or_else(|| crate::error::value(format!("{symbol:?} is not a chemical element symbol")))
}

#[pymethods]
impl PyMolecule {
    /// Builds a molecule from element symbols, coordinates and `(first, second, order)` bonds.
    ///
    /// `order` is the MDL order: 1, 2, 3, or 4 for aromatic. `formal_charges` is one integer
    /// per atom (0 for none). `properties` become SDF data fields.
    #[new]
    #[allow(clippy::needless_pass_by_value)]
    #[pyo3(signature = (elements, coordinates, bonds=None, *, name="", formal_charges=None, properties=None))]
    fn construct(
        elements: Vec<String>,
        coordinates: &Bound<'_, PyArray2<f32>>,
        bonds: Option<&Bound<'_, PyArray2<u32>>>,
        name: &str,
        formal_charges: Option<Vec<i8>>,
        properties: Option<BTreeMap<String, String>>,
    ) -> PyResult<Self> {
        let coordinates = coordinates.readonly();
        let view = coordinates.as_array();
        if view.ncols() != 3 || view.nrows() != elements.len() {
            return Err(crate::error::value(
                "coordinates must have one row of three values per element",
            ));
        }
        let atoms = elements
            .iter()
            .zip(view.rows())
            .map(|(symbol, row)| {
                Ok(MolAtom {
                    element: element_of(symbol)?,
                    position: [row[0], row[1], row[2]],
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        let bonds = match bonds {
            Some(table) => {
                let table = table.readonly();
                let view = table.as_array();
                if view.ncols() != 3 {
                    return Err(crate::error::value(
                        "bonds must have three columns: first, second, order",
                    ));
                }
                view.rows()
                    .into_iter()
                    .map(|row| {
                        Ok(MolBond {
                            first: row[0] as usize,
                            second: row[1] as usize,
                            order: u8::try_from(row[2]).map_err(|_| {
                                crate::error::value("a bond order must be 1, 2, 3 or 4")
                            })?,
                        })
                    })
                    .collect::<PyResult<Vec<_>>>()?
            }
            None => Vec::new(),
        };
        let charges = match formal_charges {
            Some(charges) if charges.len() == atoms.len() => charges,
            Some(_) => {
                return Err(crate::error::value(
                    "formal_charges must have one integer per atom",
                ));
            }
            None => vec![0; atoms.len()],
        };
        let atom_metadata = charges
            .iter()
            .map(|&charge| MolAtomMetadata {
                formal_charge: (charge != 0).then_some(charge),
                ..MolAtomMetadata::default()
            })
            .collect();
        let bond_metadata = vec![MolBondMetadata::default(); bonds.len()];
        let properties = properties
            .into_iter()
            .flatten()
            .map(|(name, value)| SdfProperty {
                name: name.into(),
                value: value.into(),
            })
            .collect();
        Ok(Self::new(MolRecord {
            name: name.into(),
            program: Box::default(),
            comment: Box::default(),
            version: MolVersion::V2000,
            molecule: Molecule { atoms, bonds },
            atom_metadata,
            bond_metadata,
            properties,
        }))
    }

    /// The first header line.
    #[getter]
    fn name(&self) -> &str {
        &self.record.name
    }

    /// The second header line: the program that wrote the record.
    #[getter]
    fn program(&self) -> &str {
        &self.record.program
    }

    /// The third header line.
    #[getter]
    fn comment(&self) -> &str {
        &self.record.comment
    }

    /// `v2000` or `v3000`.
    #[getter]
    fn version(&self) -> &'static str {
        match self.record.version {
            MolVersion::V2000 => "v2000",
            MolVersion::V3000 => "v3000",
        }
    }

    #[getter]
    fn atom_count(&self) -> usize {
        self.record.molecule.atoms.len()
    }

    #[getter]
    fn bond_count(&self) -> usize {
        self.record.molecule.bonds.len()
    }

    /// Element symbols, in file order.
    #[getter]
    fn elements(&self) -> Vec<&'static str> {
        self.record
            .molecule
            .atoms
            .iter()
            .map(|atom| atom.element.symbol())
            .collect()
    }

    /// Atom positions, `(n, 3)`.
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let flat: Vec<f32> = self
            .record
            .molecule
            .atoms
            .iter()
            .flat_map(|atom| atom.position)
            .collect();
        let array = flat
            .into_pyarray(py)
            .reshape((self.record.molecule.atoms.len(), 3))?;
        array.readwrite().make_nonwriteable();
        Ok(array)
    }

    /// Bonds as `(first, second, order)` rows, zero-based; order 4 is aromatic.
    #[getter]
    fn bonds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let mut flat = Vec::with_capacity(self.record.molecule.bonds.len() * 3);
        for bond in &self.record.molecule.bonds {
            flat.push(u32::try_from(bond.first).map_err(|_| crate::error::value("bond index"))?);
            flat.push(u32::try_from(bond.second).map_err(|_| crate::error::value("bond index"))?);
            flat.push(u32::from(bond.order));
        }
        let array = flat
            .into_pyarray(py)
            .reshape((self.record.molecule.bonds.len(), 3))?;
        array.readwrite().make_nonwriteable();
        Ok(array)
    }

    /// Known MDL formal charges; an omitted charge is neutral (zero).
    #[getter]
    fn formal_charges(&self) -> Vec<i8> {
        self.record
            .atom_metadata
            .iter()
            .map(|metadata| match metadata.formal_charge {
                Some(charge) => charge,
                None => 0,
            })
            .collect()
    }

    /// The SDF data fields, in source order.
    #[getter]
    fn properties<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let fields = PyDict::new(py);
        for property in &self.record.properties {
            fields.set_item(&*property.name, &*property.value)?;
        }
        Ok(fields)
    }

    /// The molecule as a structure with source bonds and known MDL formal charges.
    /// Explicit aromatic bonds and conservatively perceived conjugated circuits
    /// annotate aromatic atoms. Unproved aromaticity remains unknown.
    fn to_structure(&self) -> PyResult<PyStructure> {
        mol_record_to_structure(&self.record)
            .map(|structure| PyStructure::new(molframe::Structure::from(structure)))
            .map_err(crate::error::kernel)
    }

    /// The record as a MOL block.
    fn to_mol(&self) -> PyResult<String> {
        write_mol(&self.record).map_err(crate::error::kernel)
    }

    /// The record as one SDF entry, data fields included.
    fn to_sdf(&self) -> PyResult<String> {
        write_sdf(std::slice::from_ref(&*self.record)).map_err(crate::error::kernel)
    }

    fn __repr__(&self) -> String {
        format!(
            "Molecule(name={:?}, atoms={}, bonds={})",
            &*self.record.name,
            self.record.molecule.atoms.len(),
            self.record.molecule.bonds.len()
        )
    }
}

/// Every record of an SDF file, with its data fields.
#[pyfunction]
fn read_sdf(py: Python<'_>, source: &Bound<'_, PyAny>) -> PyResult<Vec<PyMolecule>> {
    let text = text_of(source)?;
    let records = py
        .detach(|| parse_sdf_records(&text))
        .map_err(crate::error::kernel)?;
    Ok(records.into_iter().map(PyMolecule::new).collect())
}

/// One MOL block.
#[pyfunction]
fn read_mol(py: Python<'_>, source: &Bound<'_, PyAny>) -> PyResult<PyMolecule> {
    let text = text_of(source)?;
    let record = py
        .detach(|| parse_mol_record(&text))
        .map_err(crate::error::kernel)?;
    Ok(PyMolecule::new(record))
}

/// SDF text for the molecules, each followed by its data fields.
#[pyfunction]
#[pyo3(name = "write_sdf")]
fn write_sdf_text(py: Python<'_>, molecules: Vec<PyRef<'_, PyMolecule>>) -> PyResult<String> {
    let records: Vec<Arc<MolRecord>> = molecules.iter().map(|m| Arc::clone(&m.record)).collect();
    drop(molecules);
    py.detach(|| {
        let owned: Vec<MolRecord> = records.iter().map(|record| (**record).clone()).collect();
        write_sdf(&owned)
    })
    .map_err(crate::error::kernel)
}

/// The graph of a structure as a molecule: elements, positions and bonds.
#[pyfunction]
#[pyo3(signature = (structure, *, name=""))]
fn molecule(py: Python<'_>, structure: &PyStructure, name: &str) -> PyResult<PyMolecule> {
    let source = structure.inner.clone();
    let graph = py
        .detach(move || structure_to_molecule(source.engine()))
        .map_err(crate::error::kernel)?;
    let atoms = graph.atoms.len();
    let bonds = graph.bonds.len();
    Ok(PyMolecule::new(MolRecord {
        name: name.into(),
        program: Box::default(),
        comment: Box::default(),
        version: MolVersion::V2000,
        molecule: graph,
        atom_metadata: vec![MolAtomMetadata::default(); atoms],
        bond_metadata: vec![MolBondMetadata::default(); bonds],
        properties: Vec::new(),
    }))
}

/// Every mapping of a SMARTS pattern onto the structure, as tuples of atom indices.
///
/// The structure needs the data the pattern's primitives use (a bond graph; aromaticity,
/// formal charges or stereochemistry where the pattern asks for them), and a missing one is
/// refused rather than guessed.
#[pyfunction]
fn smarts(py: Python<'_>, structure: &PyStructure, pattern: &str) -> PyResult<Vec<Vec<usize>>> {
    let pattern = SmartsPattern::parse(pattern).map_err(crate::error::kernel)?;
    let source = structure.inner.clone();
    let found = py
        .detach(move || pattern.find_structure_matches(source.engine()))
        .map_err(crate::error::kernel)?;
    Ok(found
        .into_iter()
        .map(|found| found.atom_indices.into_vec())
        .collect())
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMolecule>()?;
    module.add_function(wrap_pyfunction!(read_sdf, module)?)?;
    module.add_function(wrap_pyfunction!(read_mol, module)?)?;
    module.add_function(wrap_pyfunction!(write_sdf_text, module)?)?;
    module.add_function(wrap_pyfunction!(molecule, module)?)?;
    module.add_function(wrap_pyfunction!(smarts, module)?)
}
