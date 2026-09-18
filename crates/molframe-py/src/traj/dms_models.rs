//! Python value objects for the public DMS model.

use molframe::traj::{DmsBond, DmsCell, DmsFrame, DmsParticle, DmsTopology, DmsVersion};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const FINITE_VALUE: &str = "DMS numeric values must be finite";
const INVALID_BOND: &str = "DMS bonds require p0 < p1 and a finite order";

#[pyclass(name = "DmsParticle", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDmsParticle {
    pub(crate) value: DmsParticle,
}

#[pymethods]
impl PyDmsParticle {
    #[new]
    #[pyo3(signature = (*, atomic_number=None, component=None, nonbonded_type=None, mass=None, charge=None, residue_id=None, residue_name=None, chain=None, segment=None, name=None, insertion=None, formal_charge=None, occupancy=None, b_factor=None, temperature_group=None, energy_group=None, ligand_group=None, bias_group=None))]
    fn new(
        atomic_number: Option<u16>,
        component: Option<i64>,
        nonbonded_type: Option<i64>,
        mass: Option<f64>,
        charge: Option<f64>,
        residue_id: Option<i64>,
        residue_name: Option<String>,
        chain: Option<String>,
        segment: Option<String>,
        name: Option<String>,
        insertion: Option<String>,
        formal_charge: Option<f64>,
        occupancy: Option<f64>,
        b_factor: Option<f64>,
        temperature_group: Option<i64>,
        energy_group: Option<i64>,
        ligand_group: Option<i64>,
        bias_group: Option<i64>,
    ) -> PyResult<Self> {
        for value in [mass, charge, formal_charge, occupancy, b_factor]
            .into_iter()
            .flatten()
        {
            if !value.is_finite() {
                return Err(PyValueError::new_err(FINITE_VALUE));
            }
        }
        Ok(Self {
            value: DmsParticle {
                atomic_number,
                component,
                nonbonded_type,
                mass,
                charge,
                residue_id,
                residue_name: residue_name.map(Into::into),
                chain: chain.map(Into::into),
                segment: segment.map(Into::into),
                name: name.map(Into::into),
                insertion: insertion.map(Into::into),
                formal_charge,
                occupancy,
                b_factor,
                temperature_group,
                energy_group,
                ligand_group,
                bias_group,
            },
        })
    }

    #[getter]
    const fn atomic_number(&self) -> Option<u16> {
        self.value.atomic_number
    }
    #[getter]
    const fn component(&self) -> Option<i64> {
        self.value.component
    }
    #[getter]
    const fn nonbonded_type(&self) -> Option<i64> {
        self.value.nonbonded_type
    }
    #[getter]
    const fn mass(&self) -> Option<f64> {
        self.value.mass
    }
    #[getter]
    const fn charge(&self) -> Option<f64> {
        self.value.charge
    }
    #[getter]
    const fn residue_id(&self) -> Option<i64> {
        self.value.residue_id
    }
    #[getter]
    fn residue_name(&self) -> Option<&str> {
        self.value.residue_name.as_deref()
    }
    #[getter]
    fn chain(&self) -> Option<&str> {
        self.value.chain.as_deref()
    }
    #[getter]
    fn segment(&self) -> Option<&str> {
        self.value.segment.as_deref()
    }
    #[getter]
    fn name(&self) -> Option<&str> {
        self.value.name.as_deref()
    }
    #[getter]
    fn insertion(&self) -> Option<&str> {
        self.value.insertion.as_deref()
    }
    #[getter]
    const fn formal_charge(&self) -> Option<f64> {
        self.value.formal_charge
    }
    #[getter]
    const fn occupancy(&self) -> Option<f64> {
        self.value.occupancy
    }
    #[getter]
    const fn b_factor(&self) -> Option<f64> {
        self.value.b_factor
    }
    #[getter]
    const fn temperature_group(&self) -> Option<i64> {
        self.value.temperature_group
    }
    #[getter]
    const fn energy_group(&self) -> Option<i64> {
        self.value.energy_group
    }
    #[getter]
    const fn ligand_group(&self) -> Option<i64> {
        self.value.ligand_group
    }
    #[getter]
    const fn bias_group(&self) -> Option<i64> {
        self.value.bias_group
    }
}

impl From<DmsParticle> for PyDmsParticle {
    fn from(value: DmsParticle) -> Self {
        Self { value }
    }
}

impl From<PyDmsParticle> for DmsParticle {
    fn from(value: PyDmsParticle) -> Self {
        value.value
    }
}

#[pyclass(name = "DmsVersion", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDmsVersion {
    pub(crate) value: DmsVersion,
}

#[pymethods]
impl PyDmsVersion {
    #[new]
    fn new(major: u32, minor: u32) -> Self {
        Self {
            value: DmsVersion { major, minor },
        }
    }

    #[getter]
    const fn major(&self) -> u32 {
        self.value.major
    }

    #[getter]
    const fn minor(&self) -> u32 {
        self.value.minor
    }
}

impl From<DmsVersion> for PyDmsVersion {
    fn from(value: DmsVersion) -> Self {
        Self { value }
    }
}

impl From<PyDmsVersion> for DmsVersion {
    fn from(value: PyDmsVersion) -> Self {
        value.value
    }
}

#[pyclass(name = "DmsBond", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDmsBond {
    pub(crate) value: DmsBond,
}

#[pymethods]
impl PyDmsBond {
    #[new]
    fn new(p0: u32, p1: u32, order: f64) -> PyResult<Self> {
        if p0 >= p1 || !order.is_finite() {
            return Err(PyValueError::new_err(INVALID_BOND));
        }
        Ok(Self {
            value: DmsBond { p0, p1, order },
        })
    }

    #[getter]
    const fn p0(&self) -> u32 {
        self.value.p0
    }
    #[getter]
    const fn p1(&self) -> u32 {
        self.value.p1
    }
    #[getter]
    const fn order(&self) -> f64 {
        self.value.order
    }
}

impl From<DmsBond> for PyDmsBond {
    fn from(value: DmsBond) -> Self {
        Self { value }
    }
}

impl From<PyDmsBond> for DmsBond {
    fn from(value: PyDmsBond) -> Self {
        value.value
    }
}

#[pyclass(name = "DmsCell", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDmsCell {
    pub(crate) value: DmsCell,
}

#[pymethods]
impl PyDmsCell {
    #[new]
    fn new(vectors: [[f64; 3]; 3]) -> PyResult<Self> {
        if vectors.iter().flatten().any(|value| !value.is_finite()) {
            return Err(PyValueError::new_err(FINITE_VALUE));
        }
        Ok(Self {
            value: DmsCell { vectors },
        })
    }

    #[getter]
    const fn vectors(&self) -> [[f64; 3]; 3] {
        self.value.vectors
    }
}

impl From<DmsCell> for PyDmsCell {
    fn from(value: DmsCell) -> Self {
        Self { value }
    }
}

impl From<PyDmsCell> for DmsCell {
    fn from(value: PyDmsCell) -> Self {
        value.value
    }
}

#[pyclass(name = "DmsFrame", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDmsFrame {
    pub(crate) value: DmsFrame,
}

#[pymethods]
impl PyDmsFrame {
    #[new]
    #[pyo3(signature = (positions, *, velocities=None, cell=None))]
    fn new(
        positions: Vec<[f64; 3]>,
        velocities: Option<Vec<Option<[f64; 3]>>>,
        cell: Option<PyDmsCell>,
    ) -> PyResult<Self> {
        if positions.iter().flatten().any(|value| !value.is_finite()) {
            return Err(PyValueError::new_err(FINITE_VALUE));
        }
        let velocities = match velocities {
            Some(values) if values.len() == positions.len() => values,
            Some(_) => {
                return Err(PyValueError::new_err(
                    "DMS velocities must have one entry per position",
                ));
            }
            None => vec![None; positions.len()],
        };
        if velocities
            .iter()
            .flatten()
            .flatten()
            .any(|value| !value.is_finite())
        {
            return Err(PyValueError::new_err(FINITE_VALUE));
        }
        Ok(Self {
            value: DmsFrame {
                positions,
                velocities,
                cell: cell.map(Into::into),
            },
        })
    }

    #[getter]
    fn positions(&self) -> Vec<[f64; 3]> {
        self.value.positions.clone()
    }
    #[getter]
    fn velocities(&self) -> Vec<Option<[f64; 3]>> {
        self.value.velocities.clone()
    }
    #[getter]
    fn cell(&self) -> Option<PyDmsCell> {
        self.value.cell.map(Into::into)
    }
}

impl From<DmsFrame> for PyDmsFrame {
    fn from(value: DmsFrame) -> Self {
        Self { value }
    }
}

impl From<PyDmsFrame> for DmsFrame {
    fn from(value: PyDmsFrame) -> Self {
        value.value
    }
}

#[pyclass(name = "DmsTopology", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDmsTopology {
    pub(crate) value: DmsTopology,
}

#[pymethods]
impl PyDmsTopology {
    #[new]
    fn new(particles: Vec<PyDmsParticle>, bonds: Vec<PyDmsBond>) -> PyResult<Self> {
        let particles = particles
            .into_iter()
            .map(DmsParticle::from)
            .collect::<Vec<_>>();
        let bonds = bonds.into_iter().map(DmsBond::from).collect::<Vec<_>>();
        for bond in &bonds {
            if usize::try_from(bond.p1).map_or(true, |index| index >= particles.len()) {
                return Err(PyValueError::new_err(
                    "DMS bond particle identifiers must address the topology",
                ));
            }
        }
        Ok(Self {
            value: DmsTopology { particles, bonds },
        })
    }

    #[getter]
    fn particles(&self) -> Vec<PyDmsParticle> {
        self.value
            .particles
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }
    #[getter]
    fn bonds(&self) -> Vec<PyDmsBond> {
        self.value.bonds.iter().copied().map(Into::into).collect()
    }
}

impl From<DmsTopology> for PyDmsTopology {
    fn from(value: DmsTopology) -> Self {
        Self { value }
    }
}

impl From<PyDmsTopology> for DmsTopology {
    fn from(value: PyDmsTopology) -> Self {
        value.value
    }
}
