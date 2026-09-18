//! Format-specific trajectory projections.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "XyzAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyXyzAtom {
    #[pyo3(get)]
    pub(crate) element: String,
    #[pyo3(get)]
    pub(crate) position: [f32; 3],
}

#[pymethods]
impl PyXyzAtom {
    #[new]
    fn new(element: String, position: [f32; 3]) -> Self {
        Self { element, position }
    }
}

#[pyclass(name = "XyzFrame", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyXyzFrame {
    #[pyo3(get)]
    pub(crate) comment: String,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyXyzAtom>,
}

#[pymethods]
impl PyXyzFrame {
    #[new]
    fn new(comment: String, atoms: Vec<PyXyzAtom>) -> Self {
        Self { comment, atoms }
    }
}

#[pyfunction]
pub(crate) fn parse_xyz(py: Python<'_>, text: &str) -> Option<Vec<PyXyzFrame>> {
    let text = text.to_owned();
    let frames = py.detach(move || molframe::traj::format_xyz::parse_xyz(&text));
    frames.map(|frames| {
        frames
            .into_iter()
            .map(|frame| PyXyzFrame {
                comment: frame.comment,
                atoms: frame
                    .atoms
                    .into_iter()
                    .map(|atom| PyXyzAtom {
                        element: atom.element.symbol().to_owned(),
                        position: atom.position,
                    })
                    .collect(),
            })
            .collect()
    })
}

#[pyfunction]
pub(crate) fn write_xyz(py: Python<'_>, frames: Vec<PyXyzFrame>) -> PyResult<String> {
    py.detach(move || -> PyResult<String> {
        let frames = frames
            .into_iter()
            .map(|frame| {
                let atoms = frame
                    .atoms
                    .into_iter()
                    .map(|atom| {
                        let element =
                            molframe::Element::from_symbol(&atom.element).ok_or_else(|| {
                                PyValueError::new_err(format!("unknown element: {}", atom.element))
                            })?;
                        Ok(molframe::traj::format_xyz::XyzAtom {
                            element,
                            position: atom.position,
                        })
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                Ok(molframe::traj::format_xyz::XyzFrame {
                    comment: frame.comment,
                    atoms,
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(molframe::traj::format_xyz::write_xyz(&frames))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyXyzAtom>()?;
    module.add_class::<PyXyzFrame>()?;
    module.add_function(wrap_pyfunction!(parse_xyz, module)?)?;
    module.add_function(wrap_pyfunction!(write_xyz, module)?)?;
    Ok(())
}
