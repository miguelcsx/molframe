//! Density maps: reading and writing MRC/CCP4, their statistics and interpolation.

use crate::crystal::PyUnitCell;
use molframe::crystal::{DensityMap, MapBoundary};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArray3, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::path::PathBuf;
use std::sync::Arc;

fn boundary(word: &str) -> PyResult<MapBoundary> {
    match word {
        "missing" => Ok(MapBoundary::Missing),
        "periodic" => Ok(MapBoundary::Periodic),
        other => Err(crate::error::value(format!(
            "{other:?} is not a boundary; the boundaries are \"missing\" and \"periodic\""
        ))),
    }
}

/// A scalar density grid: values on a regular grid of a unit cell.
#[derive(Clone, Debug)]
#[pyclass(
    name = "DensityMap",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PyDensityMap {
    inner: Arc<DensityMap>,
}

/// Summary of a set of density values.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "MapStatistics",
    frozen,
    skip_from_py_object,
    module = "molframe.crystal"
)]
pub(crate) struct PyMapStatistics {
    inner: molframe::crystal::MapStatistics,
}

#[pymethods]
impl PyMapStatistics {
    /// Voxels included.
    #[getter]
    const fn count(&self) -> usize {
        self.inner.count
    }

    #[getter]
    const fn minimum(&self) -> f32 {
        self.inner.minimum
    }

    #[getter]
    const fn maximum(&self) -> f32 {
        self.inner.maximum
    }

    #[getter]
    const fn mean(&self) -> f64 {
        self.inner.mean
    }

    /// Population standard deviation.
    #[getter]
    const fn sigma(&self) -> f64 {
        self.inner.sigma
    }
}

#[pymethods]
impl PyDensityMap {
    /// Builds a map from `values` of shape `(nz, ny, nx)` (X varies fastest) on `cell`.
    ///
    /// `sampling` is the number of grid intervals the cell is divided into along x, y, z; it
    /// defaults to the grid's own size. `starts` is the first stored grid point, `origin` the
    /// real-space origin in ångström, `space_group` the International Tables number (0 for
    /// none).
    #[new]
    #[pyo3(signature = (
        values,
        cell,
        *,
        sampling=None,
        starts=[0, 0, 0],
        origin=[0.0, 0.0, 0.0],
        space_group=0,
        labels=Vec::new(),
    ))]
    #[allow(clippy::needless_pass_by_value)]
    fn new(
        values: &Bound<'_, PyArray3<f32>>,
        cell: &PyUnitCell,
        sampling: Option<[usize; 3]>,
        starts: [i32; 3],
        origin: [f64; 3],
        space_group: i32,
        labels: Vec<String>,
    ) -> Self {
        let readonly = values.readonly();
        let shape = readonly.shape();
        let dimensions = [shape[2], shape[1], shape[0]];
        let flat: Vec<f32> = readonly.as_array().iter().copied().collect();
        Self {
            inner: Arc::new(DensityMap {
                dimensions,
                starts,
                sampling: match sampling {
                    Some(sampling) => sampling,
                    None => dimensions,
                },
                cell: cell.parameters(),
                origin,
                space_group,
                labels: labels.into_iter().map(String::into_boxed_str).collect(),
                extended_header: Vec::new(),
                values: flat,
            }),
        }
    }

    /// Grid points stored along x, y, z.
    #[getter]
    fn dimensions(&self) -> (usize, usize, usize) {
        let [x, y, z] = self.inner.dimensions;
        (x, y, z)
    }

    /// The first stored grid point along x, y, z.
    #[getter]
    fn starts(&self) -> (i32, i32, i32) {
        let [x, y, z] = self.inner.starts;
        (x, y, z)
    }

    /// Grid intervals the unit cell is divided into along x, y, z.
    #[getter]
    fn sampling(&self) -> (usize, usize, usize) {
        let [x, y, z] = self.inner.sampling;
        (x, y, z)
    }

    /// The unit cell.
    #[getter]
    fn cell(&self) -> PyResult<PyUnitCell> {
        PyUnitCell::from_cell(&self.inner.cell)
    }

    /// The real-space origin in ångström.
    #[getter]
    fn origin(&self) -> (f64, f64, f64) {
        let [x, y, z] = self.inner.origin;
        (x, y, z)
    }

    /// The International Tables space-group number; 0 when the file gave none.
    #[getter]
    fn space_group(&self) -> i32 {
        self.inner.space_group
    }

    /// The header labels.
    #[getter]
    fn labels(&self) -> Vec<String> {
        self.inner.labels.iter().map(ToString::to_string).collect()
    }

    /// The densities, `(nz, ny, nx)`: x varies fastest.
    #[getter]
    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f32>>> {
        let [nx, ny, nz] = self.inner.dimensions;
        let array = self
            .inner
            .values
            .clone()
            .into_pyarray(py)
            .reshape((nz, ny, nx))?;
        array.readwrite().make_nonwriteable();
        Ok(array)
    }

    /// Mean, population sigma, minimum and maximum over every voxel.
    fn statistics(&self) -> PyResult<PyMapStatistics> {
        self.inner
            .statistics()
            .map(|inner| PyMapStatistics { inner })
            .map_err(crate::error::kernel)
    }

    /// The same over the voxels where `mask` (shape `(nz, ny, nx)`) is true.
    fn masked_statistics(&self, mask: &Bound<'_, PyArray3<bool>>) -> PyResult<PyMapStatistics> {
        let flat: Vec<bool> = mask.readonly().as_array().iter().copied().collect();
        self.inner
            .masked_statistics(&flat)
            .map(|inner| PyMapStatistics { inner })
            .map_err(crate::error::kernel)
    }

    /// Counts of voxels in `bins` equal intervals of `[minimum, maximum]`; values outside the
    /// range are left out.
    fn histogram<'py>(
        &self,
        py: Python<'py>,
        bins: usize,
        minimum: f32,
        maximum: f32,
    ) -> PyResult<Bound<'py, PyArray1<u64>>> {
        let found = self
            .inner
            .histogram(bins, minimum, maximum)
            .map_err(crate::error::kernel)?;
        let counts: Vec<u64> = found.counts.iter().map(|&count| count as u64).collect();
        Ok(counts.into_pyarray(py))
    }

    /// The density at Cartesian positions (ångström), interpolated; `nan` where there is none.
    ///
    /// `method` is `"linear"` or `"cubic"` (Catmull-Rom); `boundary` is `"missing"` (outside the
    /// stored grid there is no value) or `"periodic"` (the grid repeats).
    #[pyo3(signature = (positions, *, method="linear", boundary="missing"))]
    fn sample<'py>(
        &self,
        py: Python<'py>,
        positions: &Bound<'py, PyArray2<f64>>,
        method: &str,
        boundary: &str,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        let boundary = self::boundary(boundary)?;
        let cubic = match method {
            "linear" => false,
            "cubic" => true,
            other => {
                return Err(crate::error::value(format!(
                    "{other:?} is not a method; the methods are \"linear\" and \"cubic\""
                )));
            }
        };
        let readonly = positions.readonly();
        let view = readonly.as_array();
        if view.ncols() != 3 {
            return Err(crate::error::value("positions must have three columns"));
        }
        let points: Vec<[f64; 3]> = view
            .rows()
            .into_iter()
            .map(|row| [row[0], row[1], row[2]])
            .collect();
        let map = Arc::clone(&self.inner);
        let found = py.detach(move || {
            let sampler = map.sampler();
            points
                .iter()
                .map(|&point| {
                    let value = sampler.as_ref().and_then(|sampler| {
                        if cubic {
                            sampler.sample_cartesian_cubic(point, boundary)
                        } else {
                            sampler.sample_cartesian(point, boundary)
                        }
                    });
                    match value {
                        Some(value) => value,
                        None => f32::NAN,
                    }
                })
                .collect::<Vec<f32>>()
        });
        Ok(found.into_pyarray(py))
    }

    /// The map as MRC2014 bytes (mode 2, little-endian).
    fn to_mrc<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let map = Arc::clone(&self.inner);
        let bytes = py
            .detach(move || map.to_mrc_bytes())
            .map_err(crate::error::kernel)?;
        Ok(PyBytes::new(py, &bytes))
    }

    fn __repr__(&self) -> String {
        let [nx, ny, nz] = self.inner.dimensions;
        format!("DensityMap(dimensions=({nx}, {ny}, {nz}))")
    }
}

/// Reads an MRC2014 or CCP4 scalar map (complex modes are refused, not truncated).
#[pyfunction]
fn read_mrc(py: Python<'_>, source: &Bound<'_, PyAny>) -> PyResult<PyDensityMap> {
    let bytes: Vec<u8> = if let Ok(path) = source.extract::<PathBuf>() {
        std::fs::read(&path).map_err(|error| {
            crate::error::kernel(
                molframe::Diagnostic::new(molframe::Code::E7101)
                    .with_message(format!("{}: {error}", path.display())),
            )
        })?
    } else {
        source
            .extract::<pyo3::pybacked::PyBackedBytes>()
            .map_err(|_| crate::error::type_error("source must be a path or bytes"))?
            .to_vec()
    };
    let map = py
        .detach(move || DensityMap::from_mrc_bytes(&bytes))
        .map_err(crate::error::kernel)?;
    Ok(PyDensityMap {
        inner: Arc::new(map),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDensityMap>()?;
    module.add_class::<PyMapStatistics>()?;
    module.add_function(wrap_pyfunction!(read_mrc, module)?)
}
