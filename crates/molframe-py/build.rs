//! Link the Python extension independently of the consuming build frontend.

#![forbid(unsafe_code)]

fn main() {
    pyo3_build_config::add_extension_module_link_args();
}
