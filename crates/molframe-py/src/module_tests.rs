//! The extension module, built in process and held to what it publishes.
//!
//! Building it *is* the check. `namespace_registration::register` copies every
//! name a namespace advertises out of the root module with `getattr` and
//! propagates the miss, so a name that a namespace lists and nothing registers
//! fails the build. This is the only place that build runs before the wheel
//! exists: until it did, the first thing to notice was `import molframe`.
//!
//! Only that direction is asserted. Registering a name that no namespace
//! publishes is not a defect — it is how the root keeps a vocabulary larger than
//! the sum of its packages.

use pyo3::prelude::*;
use pyo3::types::PyModule;

#[test]
fn the_module_builds_with_every_name_its_namespaces_advertise() {
    Python::initialize();
    Python::attach(|py| {
        let module = match PyModule::new(py, "_native") {
            Ok(module) => module,
            Err(error) => panic!("test module should be created: {error}"),
        };
        if let Err(error) = super::native(&module) {
            panic!("the extension module should register: {error}");
        }
        for (namespace, exports) in super::namespace_registration::NAMESPACES {
            let package = match module.getattr(namespace) {
                Ok(package) => package,
                Err(error) => panic!("{namespace} should be a namespace: {error}"),
            };
            for name in *exports {
                assert!(
                    package.getattr(name).is_ok(),
                    "{namespace} advertises {name} without carrying it"
                );
            }
        }
    });
}
