//! mmCIF text for the crystallographic categories a structure carries.
//!
//! These are the inverse of the lowerings in this crate: a file read into a
//! structure keeps its assemblies, non-crystallographic operators and space
//! group as typed extensions, and a writer needs them back as categories so the
//! round trip does not silently drop them.

use crate::{AssemblySet, NcsCode, NcsSet, SymmetrySet};
use molframe_cif::quote_text;
use std::fmt::Write as _;

fn number(value: f64) -> String {
    // The shortest text that reads back to the same double.
    format!("{value}")
}

fn optional(value: Option<&str>) -> String {
    match value {
        Some(text) if !text.is_empty() => quote_text(text),
        _ => "?".to_owned(),
    }
}

/// `_pdbx_struct_assembly`, `_pdbx_struct_assembly_gen` and
/// `_pdbx_struct_oper_list` for an assembly set; empty when there is none.
#[must_use]
pub fn assemblies_cif(set: &AssemblySet) -> String {
    let mut out = String::new();
    if set.is_empty() {
        return out;
    }
    out.push_str(
        "loop_\n_pdbx_struct_assembly.id\n_pdbx_struct_assembly.details\n\
         _pdbx_struct_assembly.method_details\n_pdbx_struct_assembly.oligomeric_details\n\
         _pdbx_struct_assembly.oligomeric_count\n",
    );
    for assembly in set.assemblies() {
        let count = assembly
            .oligomeric
            .map_or_else(|| "?".to_owned(), |count| count.to_string());
        let _ = writeln!(
            out,
            "{} {} {} {} {count}",
            quote_text(&assembly.id),
            optional(assembly.details.as_deref()),
            optional(assembly.method.as_deref()),
            assembly.oligomeric.map_or("?", oligomeric_name),
        );
    }
    out.push_str("#\n");
    if set
        .assemblies()
        .any(|assembly| !assembly.generators.is_empty())
    {
        out.push_str(
            "loop_\n_pdbx_struct_assembly_gen.assembly_id\n\
             _pdbx_struct_assembly_gen.oper_expression\n\
             _pdbx_struct_assembly_gen.asym_id_list\n",
        );
        for assembly in set.assemblies() {
            for generator in &assembly.generators {
                let factors = generator.oper_expression.factors();
                let expression = match factors {
                    [only] if only.len() == 1 => only[0].to_string(),
                    _ => factors.iter().fold(String::new(), |mut text, factor| {
                        let ids: Vec<&str> = factor.iter().map(AsRef::as_ref).collect();
                        let _ = write!(text, "({})", ids.join(","));
                        text
                    }),
                };
                let asyms: Vec<&str> = generator.asym_ids.iter().map(AsRef::as_ref).collect();
                let _ = writeln!(
                    out,
                    "{} {} {}",
                    quote_text(&assembly.id),
                    quote_text(&expression),
                    quote_text(&asyms.join(",")),
                );
            }
        }
        out.push_str("#\n");
    }
    if set.operators().next().is_some() {
        out.push_str("loop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.type\n");
        matrix_header(&mut out, "pdbx_struct_oper_list");
        for operator in set.operators() {
            let _ = write!(out, "{} ? ", quote_text(&operator.id));
            matrix_row(
                &mut out,
                &operator.transform.rotation,
                &operator.transform.translation,
            );
        }
        out.push_str("#\n");
    }
    out
}

/// The dictionary's word for an oligomeric count; several readers require the
/// item to be present, so an unnamed count is written as unknown.
const fn oligomeric_name(count: u32) -> &'static str {
    match count {
        1 => "monomeric",
        2 => "dimeric",
        3 => "trimeric",
        4 => "tetrameric",
        5 => "pentameric",
        6 => "hexameric",
        7 => "heptameric",
        8 => "octameric",
        9 => "nonameric",
        10 => "decameric",
        11 => "undecameric",
        12 => "dodecameric",
        _ => "?",
    }
}

/// `_struct_ncs_oper` for an operator set; empty when there is none.
#[must_use]
pub fn ncs_cif(set: &NcsSet) -> String {
    let mut out = String::new();
    if set.is_empty() {
        return out;
    }
    out.push_str("loop_\n_struct_ncs_oper.id\n_struct_ncs_oper.code\n_struct_ncs_oper.details\n");
    matrix_header(&mut out, "struct_ncs_oper");
    for operator in set.operators() {
        let code = match operator.code {
            NcsCode::Given => "given",
            NcsCode::Generate => "generate",
        };
        let _ = write!(
            out,
            "{} {code} {} ",
            quote_text(&operator.id),
            optional(operator.details.as_deref()),
        );
        matrix_row(
            &mut out,
            &operator.transform.matrix,
            &operator.transform.translation,
        );
    }
    out.push_str("#\n");
    out
}

/// Space-group identity as `_symmetry` and `_space_group`; empty when unknown.
#[must_use]
pub fn symmetry_cif(set: &SymmetrySet) -> String {
    let mut out = String::new();
    if let Some(name) = set.hermann_mauguin.as_deref() {
        let _ = writeln!(out, "_symmetry.space_group_name_H-M   {}", quote_text(name));
    }
    if let Some(number) = set.international_number {
        let _ = writeln!(out, "_symmetry.Int_Tables_number   {number}");
    }
    if !out.is_empty() {
        out.push_str("#\n");
    }
    if let Some(hall) = set.hall.as_deref() {
        let _ = writeln!(out, "_space_group.name_Hall   {}\n#", quote_text(hall));
    }
    out
}

fn matrix_header(out: &mut String, category: &str) {
    for row in 1..=3 {
        for column in 1..=3 {
            let _ = writeln!(out, "_{category}.matrix[{row}][{column}]");
        }
    }
    for row in 1..=3 {
        let _ = writeln!(out, "_{category}.vector[{row}]");
    }
}

fn matrix_row(out: &mut String, matrix: &[[f64; 3]; 3], translation: &[f64; 3]) {
    for value in matrix.iter().flatten().chain(translation) {
        let _ = write!(out, "{} ", number(*value));
    }
    // Replace the trailing separator with the row end.
    out.pop();
    out.push('\n');
}

#[cfg(test)]
#[path = "cif_text_tests.rs"]
mod tests;
