use super::super::{chemistry, coordinates, geometry, selection, spatial, structure, surface};
use pyo3::prelude::*;
use pyo3::types::{PyAnyMethods, PyDict, PyModule};

#[test]
fn public_from_dict_rejects_foreign_operation_tags_before_reading_fields() {
    Python::initialize();
    Python::attach(|py| {
        let module = match PyModule::new(py, "declarative") {
            Ok(module) => module,
            Err(error) => panic!("test module should be created: {error}"),
        };
        register_operation_classes(&module);
        for class in [
            "Rmsd",
            "SelectQuery",
            "InferBonds",
            "Centroid",
            "NeighborPairs",
            "Sasa",
            "BasePairs",
        ] {
            reject_foreign_tag(py, &module, class);
        }
    });
}

fn register_operation_classes(module: &Bound<'_, PyModule>) {
    for registration in [
        coordinates::register as fn(&Bound<'_, PyModule>) -> PyResult<()>,
        selection::register,
        chemistry::register,
        geometry::register,
        spatial::register,
        surface::register,
        structure::register,
    ] {
        if let Err(error) = registration(module) {
            panic!("declarative operation should register: {error}");
        }
    }
}

fn reject_foreign_tag(py: Python<'_>, module: &Bound<'_, PyModule>, class: &str) {
    let configuration = PyDict::new(py);
    if let Err(error) = configuration.set_item("operation", "other_operation") {
        panic!("test configuration should accept an operation tag: {error}");
    }
    let operation = match module.getattr(class) {
        Ok(operation) => operation,
        Err(error) => panic!("{class} should be registered: {error}"),
    };
    let Err(error) = operation.call_method1("from_dict", (&configuration,)) else {
        panic!("{class}.from_dict must reject a foreign operation tag");
    };
    assert!(
        error.to_string().contains("expected serialized operation"),
        "{class}.from_dict returned an unclear error: {error}"
    );
}
