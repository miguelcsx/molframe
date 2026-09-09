//! Zero-copy Python views over the native XVG time-series parser.

use numpy::ndarray::ArrayView1;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

#[pyclass(name = "XvgSeries", frozen, skip_from_py_object)]
pub(crate) struct PyXvgSeries {
    inner: pdbiox::traj::xvg::Series,
}

#[pymethods]
impl PyXvgSeries {
    #[getter]
    fn title(&self) -> Option<&str> {
        self.inner.title.as_deref()
    }

    #[getter]
    fn abscissa_label(&self) -> Option<&str> {
        self.inner.abscissa_label.as_deref()
    }

    #[getter]
    fn ordinate_label(&self) -> Option<&str> {
        self.inner.ordinate_label.as_deref()
    }

    #[getter]
    fn legends(&self) -> Vec<&str> {
        self.inner.legends.iter().map(AsRef::as_ref).collect()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[getter]
    fn abscissa<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        readonly(slf, &slf.borrow().inner.abscissa)
    }

    fn column<'py>(slf: &Bound<'py, Self>, index: usize) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let series = slf.borrow();
        let values = series
            .inner
            .columns
            .get(index)
            .ok_or_else(|| PyIndexError::new_err("XVG column index is out of range"))?;
        Ok(readonly(slf, values))
    }

    fn column_by_legend<'py>(
        slf: &Bound<'py, Self>,
        legend: &str,
    ) -> Option<Bound<'py, PyArray1<f64>>> {
        let series = slf.borrow();
        series
            .inner
            .column_by_legend(legend)
            .map(|values| readonly(slf, values))
    }

    fn sample(&self, column: usize, abscissa: f64) -> Option<f64> {
        self.inner.sample(column, abscissa)
    }
}

#[pyfunction]
fn read_xvg(py: Python<'_>, data: &Bound<'_, PyAny>) -> PyResult<PyXvgSeries> {
    let bytes = match data.cast::<PyBytes>() {
        Ok(value) => value.as_bytes().to_vec(),
        Err(_) => data.cast::<PyString>()?.to_str()?.as_bytes().to_vec(),
    };
    py.detach(move || pdbiox::traj::xvg::read_xvg(&bytes))
        .map(|inner| PyXvgSeries { inner })
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

fn readonly<'py>(owner: &Bound<'py, PyXvgSeries>, values: &[f64]) -> Bound<'py, PyArray1<f64>> {
    let view = ArrayView1::from(values);
    // SAFETY: the NumPy base object owns `values`; PyXvgSeries is frozen, so
    // neither its vectors nor their allocations move while the view exists.
    let array = unsafe { PyArray1::borrow_from_array(&view, owner.clone().into_any()) };
    let _readonly = array.readwrite().make_nonwriteable();
    array
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyXvgSeries>()?;
    module.add_function(wrap_pyfunction!(read_xvg, module)?)
}
