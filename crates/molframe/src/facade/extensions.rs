//! Optional CIF-domain projections attached by the public facade.

#[cfg(feature = "modelcif")]
use molframe_core::Code;
use molframe_core::io::{ReadOptions, ReadResult};
use molframe_core::{Diagnostic, Structure};

pub(super) fn keep_non_model_extension_category(category: &str) -> bool {
    #[cfg(feature = "crystal")]
    if matches!(
        category,
        "pdbx_struct_oper_list"
            | "pdbx_struct_assembly"
            | "pdbx_struct_assembly_gen"
            | "struct_ncs_oper"
            | "space_group"
            | "symmetry"
            | "space_group_symop"
            | "symmetry_equiv"
    ) {
        return true;
    }
    #[cfg(feature = "chemistry")]
    if matches!(category, "chem_comp" | "chem_comp_atom" | "chem_comp_bond") {
        return true;
    }
    #[cfg(not(any(feature = "crystal", feature = "chemistry")))]
    let _ = category;
    false
}

/// Attaches extensions from a document supplied by an explicitly materialising format.
///
/// The component tables are moved out of `document` rather than copied.
pub(super) fn attach_materialized_cif_metadata(
    document: &mut molframe_cif::Document,
    structure: Structure,
    findings: Vec<Diagnostic>,
    options: &ReadOptions,
) -> ReadResult {
    #[cfg(feature = "modelcif")]
    let (structure, findings) = {
        let projected = molframe_modelcif::lower(&*document);
        attach_model(structure, findings, projected)
    };
    let (structure, findings) = attach_non_model(document, structure, findings);
    options.finish(structure, findings)
}

#[cfg(feature = "modelcif")]
pub(super) fn attach_projected_metadata(
    document: &mut molframe_cif::Document,
    structure: Structure,
    findings: Vec<Diagnostic>,
    projected: Result<
        (molframe_modelcif::ModelCif, Vec<Diagnostic>),
        molframe_modelcif::ModelCifError,
    >,
    options: &ReadOptions,
) -> ReadResult {
    let (structure, findings) = attach_model(structure, findings, projected);
    let (structure, findings) = attach_non_model(document, structure, findings);
    options.finish(structure, findings)
}

#[cfg(feature = "modelcif")]
fn attach_model(
    mut structure: Structure,
    mut findings: Vec<Diagnostic>,
    projected: Result<
        (molframe_modelcif::ModelCif, Vec<Diagnostic>),
        molframe_modelcif::ModelCifError,
    >,
) -> (Structure, Vec<Diagnostic>) {
    match projected {
        Ok((model, model_findings)) => {
            findings.extend(model_findings);
            if !model.is_empty() {
                match model.plddt_annotation(&structure) {
                    Ok(Some(column)) => {
                        let mut data = structure.data().clone();
                        let _ = data.annotations.insert(
                            molframe_core::annotation::PLDDT_ANNOTATION,
                            molframe_core::annotation::AtomAnnotation::Real(column),
                        );
                        structure = Structure::new(data);
                    }
                    Ok(None) => {}
                    Err(error) => findings.push(
                        Diagnostic::new(Code::E1901)
                            .with_message(error.to_string())
                            .in_category("ma_qa_metric_local"),
                    ),
                }
                structure = structure.with_extension(molframe_modelcif::MODEL_CIF_EXTENSION, model);
            }
        }
        Err(error) => findings.push(model_error(&error)),
    }
    (structure, findings)
}

#[cfg(feature = "modelcif")]
pub(super) fn model_error(error: &molframe_modelcif::ModelCifError) -> Diagnostic {
    Diagnostic::new(Code::E1901)
        .with_message(error.to_string())
        .in_category("ma_*")
}

fn attach_non_model(
    document: &mut molframe_cif::Document,
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    let (structure, findings) = attach_embedded_components(document, structure, findings);
    attach_crystal(&*document, structure, findings)
}

/// Keeps the component definitions the file carries, so annotation resolves
/// them before any dictionary without the caller asking.
fn attach_embedded_components(
    document: &mut molframe_cif::Document,
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    #[cfg(feature = "chemistry")]
    {
        let mut structure = structure;
        let mut findings = findings;
        match molframe_chem::take_embedded_components(document) {
            Ok(Some(provider)) => {
                structure = structure
                    .with_extension(molframe_chem::EMBEDDED_COMPONENTS_EXTENSION, provider);
            }
            Ok(None) => {}
            Err(values) => findings.push(
                Diagnostic::new(molframe_core::Code::W3206)
                    .in_category("chem_comp_atom")
                    .with_context("problems", values.len().to_string())
                    .with_context(
                        "first",
                        match values.first() {
                            Some(first) => first.to_string(),
                            None => String::new(),
                        },
                    ),
            ),
        }
        (structure, findings)
    }
    #[cfg(not(feature = "chemistry"))]
    {
        let _ = document;
        (structure, findings)
    }
}

fn attach_crystal(
    document: &molframe_cif::Document,
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    #[cfg(feature = "crystal")]
    {
        let mut structure = structure;
        let mut findings = findings;
        match molframe_xtal::lower_assemblies(document) {
            Ok(assemblies) if !assemblies.is_empty() => {
                structure =
                    structure.with_extension(molframe_xtal::ASSEMBLIES_EXTENSION, assemblies);
            }
            Ok(_) => {}
            Err(values) => findings.extend(values),
        }
        match molframe_xtal::lower_ncs(document) {
            Ok(ncs) if !ncs.is_empty() => {
                structure = structure.with_extension(molframe_xtal::NCS_EXTENSION, ncs);
            }
            Ok(_) => {}
            Err(values) => findings.extend(values),
        }
        match molframe_xtal::lower_symmetry(document) {
            Ok(symmetry) if !symmetry.is_empty() => {
                structure = structure.with_extension(molframe_xtal::SYMMETRY_EXTENSION, symmetry);
            }
            Ok(_) => {}
            Err(values) => findings.extend(values),
        }
        (structure, findings)
    }
    #[cfg(not(feature = "crystal"))]
    {
        let _ = document;
        (structure, findings)
    }
}
