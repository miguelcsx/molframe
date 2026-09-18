//! Native frame transforms composed by [`PipelineReader`].

use super::generic::PyMemoryReader;
use super::reader_types::PyTimestep;
use numpy::{PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[derive(Clone, Debug)]
enum TransformSpec {
    Rigid(molframe::traj::RigidTransform),
    Center(Box<[usize]>, [f64; 3]),
    Fit(Box<[usize]>, Box<[[f32; 3]]>),
    Wrap(Option<Vec<Box<[usize]>>>),
    Unwrap(Box<[(usize, usize)]>),
}

#[pyclass(name = "RigidTransform", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyRigidTransform {
    inner: molframe::traj::RigidTransform,
}

#[pymethods]
impl PyRigidTransform {
    #[new]
    fn new(offset: [f64; 3]) -> Self {
        Self {
            inner: molframe::traj::RigidTransform::translation(offset),
        }
    }

    #[staticmethod]
    fn translation(offset: [f64; 3]) -> Self {
        Self::new(offset)
    }

    #[getter]
    const fn name(&self) -> &'static str {
        "rigid"
    }
}

#[pyclass(name = "Center", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCenter {
    atoms: Box<[usize]>,
    target: [f64; 3],
}

#[pymethods]
impl PyCenter {
    #[staticmethod]
    fn geometric(atoms: Vec<usize>, target: [f64; 3]) -> Self {
        Self {
            atoms: atoms.into_boxed_slice(),
            target,
        }
    }

    #[getter]
    const fn name(&self) -> &'static str {
        "center_geometric"
    }
}

#[pyclass(name = "Fit", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFit {
    atoms: Box<[usize]>,
    reference: Box<[[f32; 3]]>,
}

#[pymethods]
impl PyFit {
    #[new]
    fn new(atoms: Vec<usize>, reference: PyReadonlyArray2<'_, f32>) -> PyResult<Self> {
        let reference = triples(reference)?;
        Ok(Self {
            atoms: atoms.into_boxed_slice(),
            reference,
        })
    }

    #[getter]
    const fn name(&self) -> &'static str {
        "fit"
    }
}

#[pyclass(name = "Wrap", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyWrap {
    groups: Option<Vec<Box<[usize]>>>,
}

#[pymethods]
impl PyWrap {
    #[staticmethod]
    fn atoms() -> Self {
        Self { groups: None }
    }

    #[staticmethod]
    fn groups(groups: Vec<Vec<usize>>) -> Self {
        Self {
            groups: Some(groups.into_iter().map(Vec::into_boxed_slice).collect()),
        }
    }
}

#[pyclass(name = "Unwrap", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyUnwrap {
    bonds: Box<[(usize, usize)]>,
}

#[pymethods]
impl PyUnwrap {
    #[staticmethod]
    fn molecules(bonds: Vec<(usize, usize)>) -> Self {
        Self {
            bonds: bonds.into_boxed_slice(),
        }
    }
}

#[pyclass(name = "PipelineReader", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPipelineReader {
    source: PyMemoryReader,
    transforms: Vec<TransformSpec>,
    cursor: usize,
}

#[pymethods]
impl PyPipelineReader {
    #[new]
    fn new(source: PyMemoryReader) -> Self {
        Self {
            source,
            transforms: Vec::new(),
            cursor: 0,
        }
    }

    fn then(&mut self, transform: &Bound<'_, PyAny>) -> PyResult<()> {
        self.transforms.push(transform_spec(transform)?);
        Ok(())
    }

    #[getter]
    fn transform_names(&self) -> Vec<&'static str> {
        self.transforms.iter().map(transform_name).collect()
    }

    #[getter]
    const fn n_atoms(&self) -> usize {
        self.source.atom_count
    }

    fn read_next(&mut self) -> PyResult<Option<PyTimestep>> {
        let Some(mut frame) = self.source.read_next_frame() else {
            return Ok(None);
        };
        frame.frame = self.cursor;
        self.cursor += 1;
        apply_transforms(&mut frame, &self.transforms)?;
        Ok(Some(frame))
    }

    fn read_all(&mut self) -> PyResult<Vec<PyTimestep>> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_next()? {
            frames.push(frame);
        }
        Ok(frames)
    }

    #[pyo3(name = "into_inner")]
    fn take_reader(&mut self) -> PyMemoryReader {
        std::mem::replace(
            &mut self.source,
            PyMemoryReader {
                frames: Vec::new(),
                cursor: 0,
                atom_count: 0,
            },
        )
    }
}

fn triples(values: PyReadonlyArray2<'_, f32>) -> PyResult<Box<[[f32; 3]]>> {
    if values.shape().get(1) != Some(&3) {
        return Err(PyValueError::new_err("coordinates must have shape (n, 3)"));
    }
    let result = values
        .as_array()
        .outer_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    drop(values);
    Ok(result)
}

fn transform_spec(value: &Bound<'_, PyAny>) -> PyResult<TransformSpec> {
    if value.is_instance_of::<PyRigidTransform>() {
        return Ok(TransformSpec::Rigid(
            value.extract::<PyRef<'_, PyRigidTransform>>()?.inner,
        ));
    }
    if value.is_instance_of::<PyCenter>() {
        let value = value.extract::<PyRef<'_, PyCenter>>()?;
        return Ok(TransformSpec::Center(value.atoms.clone(), value.target));
    }
    if value.is_instance_of::<PyFit>() {
        let value = value.extract::<PyRef<'_, PyFit>>()?;
        return Ok(TransformSpec::Fit(
            value.atoms.clone(),
            value.reference.clone(),
        ));
    }
    if value.is_instance_of::<PyWrap>() {
        let value = value.extract::<PyRef<'_, PyWrap>>()?;
        return Ok(TransformSpec::Wrap(value.groups.clone()));
    }
    if value.is_instance_of::<PyUnwrap>() {
        let value = value.extract::<PyRef<'_, PyUnwrap>>()?;
        return Ok(TransformSpec::Unwrap(value.bonds.clone()));
    }
    Err(PyValueError::new_err(
        "pipeline transform must be RigidTransform, Center, Fit, Wrap or Unwrap",
    ))
}

fn transform_name(value: &TransformSpec) -> &'static str {
    match value {
        TransformSpec::Rigid(_) => "rigid",
        TransformSpec::Center(_, _) => "center_geometric",
        TransformSpec::Fit(_, _) => "fit",
        TransformSpec::Wrap(Some(_)) => "wrap_groups",
        TransformSpec::Wrap(None) => "wrap_atoms",
        TransformSpec::Unwrap(_) => "make_molecules_whole",
    }
}

fn apply_transforms(frame: &mut PyTimestep, transforms: &[TransformSpec]) -> PyResult<()> {
    let mut native: molframe::traj::Timestep = frame.clone().try_into()?;
    for transform in transforms {
        use molframe::traj::FrameTransform;
        match transform {
            TransformSpec::Rigid(value) => value.apply(&mut native),
            TransformSpec::Center(atoms, target) => {
                molframe::traj::Center::geometric(atoms.clone(), *target).apply(&mut native)
            }
            TransformSpec::Fit(atoms, reference) => {
                molframe::traj::Fit::new(atoms.clone(), reference.clone()).apply(&mut native)
            }
            TransformSpec::Wrap(groups) => match groups {
                Some(groups) => molframe::traj::Wrap::groups(groups.clone()).apply(&mut native),
                None => molframe::traj::Wrap::atoms().apply(&mut native),
            },
            TransformSpec::Unwrap(bonds) => {
                molframe::traj::Unwrap::molecules(bonds.clone()).apply(&mut native)
            }
        }
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    }
    *frame = native.try_into()?;
    Ok(())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyRigidTransform>()?;
    module.add_class::<PyCenter>()?;
    module.add_class::<PyFit>()?;
    module.add_class::<PyWrap>()?;
    module.add_class::<PyUnwrap>()?;
    module.add_class::<PyPipelineReader>()?;
    Ok(())
}
