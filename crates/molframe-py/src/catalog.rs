//! Single checked-in operation catalog used by registration and parity checks.

use pyo3::prelude::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Capability {
    pub(crate) name: &'static str,
    pub(crate) domain: &'static str,
    pub(crate) feature: &'static str,
    pub(crate) inputs: &'static str,
    pub(crate) result: &'static str,
    pub(crate) policy: bool,
    pub(crate) eager: bool,
    pub(crate) workflow: bool,
    pub(crate) execution_needs: &'static str,
    pub(crate) cost: molframe::Cost,
}

pub(crate) const CAPABILITIES: &[Capability] = &[
    Capability {
        name: "centroid",
        domain: "geometry",
        feature: "geometry",
        inputs: "coordinates",
        result: "coordinate",
        policy: false,
        eager: true,
        workflow: true,
        execution_needs: "cpu",
        cost: molframe::Cost::Borrow,
    },
    Capability {
        name: "distance_matrix",
        domain: "geometry",
        feature: "geometry",
        inputs: "coordinates",
        result: "matrix",
        policy: false,
        eager: true,
        workflow: true,
        execution_needs: "cpu,memory",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "rmsd",
        domain: "geometry",
        feature: "geometry",
        inputs: "coordinates, coordinates",
        result: "float",
        policy: false,
        eager: true,
        workflow: true,
        execution_needs: "cpu",
        cost: molframe::Cost::Borrow,
    },
    Capability {
        name: "atom_contacts",
        domain: "analysis",
        feature: "analysis",
        inputs: "structure|selection",
        result: "ContactTable",
        policy: true,
        eager: true,
        workflow: true,
        execution_needs: "cpu,spatial,memory",
        cost: molframe::Cost::Materialize,
    },
];

pub(crate) fn validate_registration(module: &Bound<'_, PyModule>) -> PyResult<()> {
    for capability in CAPABILITIES {
        debug_assert!(!capability.inputs.is_empty());
        debug_assert!(!capability.result.is_empty());
        debug_assert!(!capability.feature.is_empty());
        debug_assert!(!capability.execution_needs.is_empty());
        debug_assert!(capability.eager || capability.workflow);
        let domain = module.getattr(capability.domain)?;
        if capability.eager {
            domain.getattr(capability.name)?;
        }
        let _ = (capability.policy, capability.cost);
    }
    Ok(())
}
