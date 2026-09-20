//! In-process registration and catalog parity checks.

use pyo3::prelude::*;
use pyo3::types::PyModule;

#[test]
fn curated_module_and_catalog_are_bidirectionally_consistent() {
    Python::initialize();
    Python::attach(|py| {
        let module = match PyModule::new(py, "_native") {
            Ok(module) => module,
            Err(error) => panic!("test module should be created: {error}"),
        };
        if let Err(error) = super::native(&module) {
            panic!("the extension module should register: {error}");
        }
        for capability in super::catalog::CAPABILITIES {
            let domain = match module.getattr(capability.domain) {
                Ok(value) => value,
                Err(error) => panic!("{} should be registered: {error}", capability.domain),
            };
            if capability.eager {
                assert!(
                    domain.getattr(capability.name).is_ok(),
                    "{}.{} is catalogued but absent",
                    capability.domain,
                    capability.name
                );
            }
        }
        for removed in ["Plan", "Batch", "AtomIndex", "StructureData", "SpatialPlan"] {
            assert!(
                module.getattr(removed).is_err(),
                "{removed} leaked into root"
            );
        }
    });
}
