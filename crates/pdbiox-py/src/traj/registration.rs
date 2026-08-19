//! Registration of trajectory-specific native entry points.

use super::formats;
use pyo3::prelude::*;

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    formats::register(module)?;
    super::data::register(module)?;
    super::generic::register(module)?;
    super::transforms::register(module)?;
    module.add_class::<super::io::PyAmberAsciiReadOptions>()?;
    module.add_class::<super::io::PyTrajectoryReadOptions>()?;
    module.add_function(wrap_pyfunction!(super::io::read_trajectory, module)?)?;
    module.add_function(wrap_pyfunction!(super::io::write_trajectory, module)?)?;
    super::kernels::register(module)?;
    super::operations::register(module)?;
    super::neighbors::register(module)?;
    super::selection::register(module)?;
    super::amber_formats::register(module)?;
    super::amber_netcdf_formats::register(module)?;
    super::binary_formats::register(module)?;
    super::container_formats::register(module)?;
    super::dlpoly_formats::register(module)?;
    super::reader_types::register(module)?;
    super::text_formats::register(module)?;
    super::topology_formats::register(module)?;
    super::gromacs_formats::register(module)?;
    super::namd_formats::register(module)?;
    super::h5md_formats::register(module)?;
    super::imd::register(module)?;
    super::text_extended::register(module)?;
    super::topology_extended::register(module)?;
    super::gamess::register(module)
}
