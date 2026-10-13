//! Metadata categories written after the entry, cell and entity blocks.
//!
//! Resolution, space group, secondary structure and caller-rendered categories
//! live here so the text and `BinaryCIF` writers emit the same set.

use super::options::CifWriteOptions;
use super::value::quoted;
use molframe_core::structure::Structure;
use std::fmt;

/// Renders every supplementary category as mmCIF text without a block header.
#[must_use]
pub fn render(structure: &Structure, options: &CifWriteOptions) -> String {
    let mut out = String::new();
    write(&mut out, structure, options);
    out
}

pub(super) fn write(out: &mut impl fmt::Write, structure: &Structure, options: &CifWriteOptions) {
    let entry = &structure.data().entry;
    let extra = match options.extra_categories() {
        Some(extra) => extra,
        None => "",
    };
    if let Some(group) = entry.space_group.as_deref()
        && !extra.contains("_symmetry.space_group_name_H-M")
        && !extra.contains("_space_group.name_H-M_alt")
    {
        let _ = writeln!(out, "_symmetry.space_group_name_H-M   {}\n#", quoted(group));
    }
    if let Some(resolution) = entry.resolution {
        // The refinement identifier is what makes a reader treat the row as a
        // refinement at all; the experiment decides which one it was.
        let method = match entry.method.as_deref() {
            Some(method) => method,
            None => "",
        };
        let refine = if method.to_ascii_uppercase().contains("NEUTRON") {
            "NEUTRON DIFFRACTION"
        } else {
            "X-RAY DIFFRACTION"
        };
        let _ = writeln!(
            out,
            "_refine.pdbx_refine_id   '{refine}'\n_refine.ls_d_res_high   {resolution:.4}\n#"
        );
    }
    super::secondary::write(out, structure);
    if !extra.is_empty() {
        let _ = out.write_str(extra);
        if !extra.ends_with('\n') {
            let _ = out.write_char('\n');
        }
    }
}
