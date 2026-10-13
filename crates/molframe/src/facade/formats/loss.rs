//! What a write to another format rewrites or cannot carry, as findings.
//!
//! A format keeps only some of what a structure holds. The writers that cannot
//! fail on that return the loss beside the output instead of dropping it, so a
//! conversion is never silent about a rewritten identifier or a discarded
//! extension.

use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::Structure;

/// Chains and residues a CIF-family write renames to keep identities distinct.
#[cfg(feature = "mmcif")]
pub(crate) fn label_rewrites(structure: &Structure) -> Vec<Diagnostic> {
    let rewrites = structure.label_rewrites();
    let mut findings = Vec::new();
    if let Some((before, after)) = rewrites.chains.first() {
        findings.push(Diagnostic::new(Code::W4102).with_message(format!(
            "{} chain(s) share a normalised label and are renamed to keep chains \
                 distinct (for example `{before}` becomes `{after}`); the depositor \
                 identifiers are unchanged",
            rewrites.chains.len()
        )));
    }
    if rewrites.residues > 0 {
        findings.push(Diagnostic::new(Code::W4102).with_message(format!(
            "{} residue(s) take a shifted normalised sequence number so insertion codes \
             stay distinct; the depositor numbers are unchanged",
            rewrites.residues
        )));
    }
    findings
}

/// The PAE matrix stored under an approximate policy, with its worst error.
///
/// Written cells are the stored approximations, so they differ from the
/// values the source file carried by up to the reported amount.
#[cfg(feature = "modelcif")]
pub(crate) fn quantised_confidence(structure: &Structure) -> Option<Diagnostic> {
    let model = structure
        .extensions()
        .get::<molframe_modelcif::ModelCif>(molframe_modelcif::MODEL_CIF_EXTENSION)?;
    let matrix = model.dense_pae()?.matrix();
    (!matrix.is_exact()).then(|| {
        Diagnostic::new(Code::W4101).with_message(format!(
            "the PAE matrix is stored quantised; written values differ from the source \
             by up to {}",
            matrix.max_quantisation_error()
        ))
    })
}

/// One finding naming every extension the target format has no place for.
pub(crate) fn uncarried_extensions(
    structure: &Structure,
    carried: &[&str],
    target: &str,
) -> Option<Diagnostic> {
    let dropped: Vec<&str> = structure
        .extensions()
        .keys()
        .filter(|key| !carried.contains(key))
        .collect();
    (!dropped.is_empty()).then(|| {
        Diagnostic::new(Code::W4101).with_message(format!(
            "extension(s) {} have no {target} counterpart and are not written",
            dropped
                .iter()
                .map(|key| format!("`{key}`"))
                .collect::<Vec<_>>()
                .join(", ")
        ))
    })
}

/// The extensions the CIF family writes back out.
#[cfg(feature = "mmcif")]
pub(crate) fn cif_carried() -> Vec<&'static str> {
    #[allow(unused_mut)]
    let mut carried = vec![molframe_core::structure::SEQUENCE_REFERENCES_EXTENSION];
    #[cfg(feature = "crystal")]
    carried.extend([
        molframe_xtal::SYMMETRY_EXTENSION,
        molframe_xtal::ASSEMBLIES_EXTENSION,
        molframe_xtal::NCS_EXTENSION,
    ]);
    #[cfg(feature = "modelcif")]
    carried.push(molframe_modelcif::MODEL_CIF_EXTENSION);
    carried
}

/// The extensions a PDB write carries; assemblies ride on the retained header
/// records, so they count only when those are present.
#[cfg(feature = "pdb")]
pub(crate) fn pdb_carried(structure: &Structure) -> Vec<&'static str> {
    #[allow(unused_mut)]
    let mut carried = vec![molframe_pdb::PDB_HEADERS_EXTENSION];
    #[cfg(feature = "crystal")]
    {
        carried.extend([
            molframe_xtal::SYMMETRY_EXTENSION,
            molframe_xtal::NCS_EXTENSION,
        ]);
        if structure
            .extensions()
            .keys()
            .any(|key| key == molframe_pdb::PDB_HEADERS_EXTENSION)
        {
            carried.push(molframe_xtal::ASSEMBLIES_EXTENSION);
        }
    }
    #[cfg(not(feature = "crystal"))]
    let _ = structure;
    carried
}
