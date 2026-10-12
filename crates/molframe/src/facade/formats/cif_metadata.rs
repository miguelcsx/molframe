//! The mmCIF write options a structure's own metadata implies.
//!
//! A structure carries assemblies, non-crystallographic operators and its space
//! group as typed extensions that the CIF crate cannot see, and it may carry no
//! data-block name at all. The facade completes the caller's options from both,
//! so a default write keeps what the file said instead of dropping it silently.

use molframe_core::structure::Structure;

/// The block name used when neither the caller nor the entry supplies one.
const FALLBACK_BLOCK_ID: &str = "molframe";

/// Options extended with the structure's own categories and a block name.
///
/// Categories and block id the caller already chose are never replaced.
pub(crate) fn complete(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> molframe_cif::CifWriteOptions {
    let mut completed = options.clone();
    if completed.extra_categories().is_none() {
        let text = extension_categories(structure);
        completed = completed.with_extra_categories(text);
    }
    if completed.block_id().is_none() && structure.data().entry.id.is_none() {
        completed = completed.with_block_id(FALLBACK_BLOCK_ID);
    }
    completed
}

#[cfg(feature = "crystal")]
fn extension_categories(structure: &Structure) -> String {
    let extensions = structure.extensions();
    let mut text = String::new();
    if let Some(set) =
        extensions.get::<molframe_xtal::SymmetrySet>(molframe_xtal::SYMMETRY_EXTENSION)
    {
        text.push_str(&molframe_xtal::symmetry_cif(set));
    }
    if let Some(set) =
        extensions.get::<molframe_xtal::AssemblySet>(molframe_xtal::ASSEMBLIES_EXTENSION)
    {
        text.push_str(&molframe_xtal::assemblies_cif(set));
    }
    if let Some(set) = extensions.get::<molframe_xtal::NcsSet>(molframe_xtal::NCS_EXTENSION) {
        text.push_str(&molframe_xtal::ncs_cif(set));
    }
    text
}

#[cfg(not(feature = "crystal"))]
fn extension_categories(_structure: &Structure) -> String {
    String::new()
}
