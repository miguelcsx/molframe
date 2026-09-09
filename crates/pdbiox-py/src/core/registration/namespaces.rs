//! Native namespace modules mirroring the facade's public crate aliases.

use pyo3::prelude::*;
use pyo3::types::PyTuple;

#[path = "namespace_analysis.rs"]
mod namespace_analysis;
#[path = "namespace_catalog.rs"]
mod namespace_catalog;
#[path = "namespace_compare.rs"]
mod namespace_compare;
#[path = "namespace_core.rs"]
mod namespace_core;
#[path = "namespace_geometry.rs"]
mod namespace_geometry;
#[path = "namespace_overlays.rs"]
mod namespace_overlays;
#[path = "namespace_platform.rs"]
mod namespace_platform;
#[path = "namespace_sequence.rs"]
mod namespace_sequence;
#[path = "trajectory_exports.rs"]
mod trajectory_exports;

use namespace_catalog::NAMESPACES;
use namespace_compare::WORKFLOW as COMPARE_WORKFLOW;
use namespace_compare::{REGION as COMPARE_REGION, SUPERPOSED as COMPARE_SUPERPOSED};
use namespace_platform::ML_ALIASES;

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    for (name, exports) in NAMESPACES {
        let namespace = PyModule::new(py, &qualified_name(module, name)?)?;
        add_exports(module, &namespace, exports)?;
        let overlays = namespace_overlays::for_namespace(name);
        add_exports(module, &namespace, overlays)?;
        for (alias, source) in namespace_overlays::aliases(name) {
            namespace.setattr(alias, module.getattr(source)?)?;
        }
        if *name == "ml" {
            for (alias, source) in ML_ALIASES {
                namespace.setattr(alias, namespace.getattr(source)?)?;
            }
        }
        if *name == "seq" || *name == "compare" {
            namespace.setattr("Column", module.getattr("AlignmentColumn")?)?;
        }
        if *name == "compare" {
            register_compare_submodules(module, &namespace)?;
        }
        if *name == "fx" {
            namespace.setattr("Evaluation", module.getattr("FxEvaluation")?)?;
        }
        if *name == "analysis" {
            namespace.setattr("DensityMap", module.getattr("_AnalysisDensityMap")?)?;
        }
        set_public_exports(
            &namespace,
            exports,
            overlays,
            namespace_overlays::alias_names(name),
            namespace_overlays::children(name),
        )?;
        add_submodule(module, &namespace)?;
    }
    module
        .getattr("seq")?
        .setattr("Code", module.getattr("SequenceCode")?)?;
    let trajectory = module.getattr("traj")?;
    let format_xyz = PyModule::new(py, "format_xyz")?;
    add_exports(
        module,
        &format_xyz,
        &["XyzAtom", "XyzFrame", "parse_xyz", "write_xyz"],
    )?;
    trajectory.setattr("format_xyz", &format_xyz)?;
    let bcif = module.getattr("bcif")?;
    bcif.setattr("read", module.getattr("read_bcif")?)?;
    bcif.setattr("read_document", module.getattr("read_bcif_document")?)?;
    bcif.setattr(
        "read_with_document",
        module.getattr("read_bcif_with_document")?,
    )?;
    bcif.setattr("write_document", module.getattr("write_bcif_document")?)?;
    Ok(())
}

fn register_compare_submodules(
    source: &Bound<'_, PyModule>,
    compare: &Bound<'_, PyModule>,
) -> PyResult<()> {
    const SUBMODULES: &[(&str, &[&str])] = &[
        ("region", COMPARE_REGION),
        ("workflow", COMPARE_WORKFLOW),
        ("superposed", COMPARE_SUPERPOSED),
    ];
    for (name, exports) in SUBMODULES {
        let submodule = PyModule::new(source.py(), &qualified_name(compare, name)?)?;
        add_exports(source, &submodule, exports)?;
        set_public_exports(&submodule, exports, &[], &[], &[])?;
        add_submodule(compare, &submodule)?;
        let root_name = format!("compare_{name}");
        source.add(&root_name, &submodule)?;
    }
    Ok(())
}

fn qualified_name(parent: &Bound<'_, PyModule>, child: &str) -> PyResult<String> {
    Ok(format!("{}.{}", parent.name()?, child))
}

fn add_submodule(parent: &Bound<'_, PyModule>, child: &Bound<'_, PyModule>) -> PyResult<()> {
    parent.add_submodule(child)?;
    let modules = PyModule::import(parent.py(), "sys")?.getattr("modules")?;
    modules.set_item(child.name()?, child)
}

pub(super) fn register_late(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.setattr("AtomSelection", module.getattr("Selection")?)?;
    let bcif = module.getattr("bcif")?;
    bcif.setattr("write_structure", module.getattr("write_bcif")?)?;
    bcif.setattr(
        "write_structure_with_options",
        module.getattr("write_bcif_with_options")?,
    )?;
    Ok(())
}

fn add_exports(
    source: &Bound<'_, PyModule>,
    target: &Bound<'_, PyModule>,
    exports: &[&str],
) -> PyResult<()> {
    for name in exports {
        let value = source.getattr(name)?;
        target.setattr(name, value)?;
    }
    Ok(())
}

fn set_public_exports(
    namespace: &Bound<'_, PyModule>,
    exports: &[&str],
    overlays: &[&str],
    aliases: &[&str],
    children: &[&str],
) -> PyResult<()> {
    let mut names =
        Vec::with_capacity(exports.len() + overlays.len() + aliases.len() + children.len());
    names.extend_from_slice(exports);
    names.extend_from_slice(overlays);
    names.extend_from_slice(aliases);
    names.extend_from_slice(children);
    let names = PyTuple::new(namespace.py(), names)?;
    namespace.setattr("__all__", names)
}
