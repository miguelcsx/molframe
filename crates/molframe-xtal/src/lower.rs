//! Lowering of PDBx/mmCIF assembly categories.

use crate::category_transform::read_affine;
use crate::{AssemblyDef, AssemblySet, Generator, OperExpression, Operator};
use molframe_cif::{Category, Document};
use molframe_core::{Code, Diagnostic};
use molframe_geom::Rigid;
use std::collections::BTreeMap;

/// Interprets the assembly categories of the first data block.
///
/// # Errors
///
/// Returns deterministic diagnostics for missing required values, duplicate
/// identifiers, malformed operator matrices, or generators naming no assembly.
pub fn lower_assemblies(document: &Document) -> Result<AssemblySet, Vec<Diagnostic>> {
    let Some(block) = document.first_block() else {
        return Ok(AssemblySet::default());
    };
    let operators = match block.category("pdbx_struct_oper_list") {
        Some(category) => lower_operators(category)?,
        None => BTreeMap::new(),
    };
    let mut assemblies = match block.category("pdbx_struct_assembly") {
        Some(category) => lower_definitions(category)?,
        None => BTreeMap::new(),
    };
    if let Some(category) = block.category("pdbx_struct_assembly_gen") {
        lower_generators(category, &mut assemblies)?;
    }
    Ok(AssemblySet::from_parts(assemblies, operators))
}

fn lower_operators(category: &Category) -> Result<BTreeMap<Box<str>, Operator>, Vec<Diagnostic>> {
    let mut operators = BTreeMap::new();
    let mut findings = Vec::new();
    for row in 0..category.row_count() {
        let Some(id) = required_identifier(category, "id", row, &mut findings) else {
            continue;
        };
        let Some(transform) = operator_transform(category, row, &mut findings) else {
            continue;
        };
        let operator = Operator {
            id: id.clone(),
            transform,
        };
        if operators.insert(id.clone(), operator).is_some() {
            findings.push(malformed(category, row, "id").with_context("value", id));
        }
    }
    finish(operators, findings)
}

fn lower_definitions(
    category: &Category,
) -> Result<BTreeMap<Box<str>, AssemblyDef>, Vec<Diagnostic>> {
    let mut assemblies = BTreeMap::new();
    let mut findings = Vec::new();
    for row in 0..category.row_count() {
        let Some(id) = required_identifier(category, "id", row, &mut findings) else {
            continue;
        };
        let oligomeric = optional_u32(category, "oligomeric_count", row, &mut findings);
        let definition = AssemblyDef {
            id: id.clone(),
            details: optional_text(category, "details", row),
            method: optional_text(category, "method_details", row),
            oligomeric,
            generators: Vec::new(),
        };
        if assemblies.insert(id.clone(), definition).is_some() {
            findings.push(malformed(category, row, "id").with_context("value", id));
        }
    }
    finish(assemblies, findings)
}

fn lower_generators(
    category: &Category,
    assemblies: &mut BTreeMap<Box<str>, AssemblyDef>,
) -> Result<(), Vec<Diagnostic>> {
    let mut findings = Vec::new();
    for row in 0..category.row_count() {
        let assembly_id = required_identifier(category, "assembly_id", row, &mut findings);
        let expression = required_identifier(category, "oper_expression", row, &mut findings);
        let asym_ids = required_identifier(category, "asym_id_list", row, &mut findings);
        let (Some(assembly_id), Some(expression), Some(asym_ids)) =
            (assembly_id, expression, asym_ids)
        else {
            continue;
        };
        let oper_expression = match OperExpression::parse(&expression) {
            Ok(expression) => expression,
            Err(finding) => {
                findings.push(finding.with_context("row", row.to_string()));
                continue;
            }
        };
        let asym_ids = parse_asym_ids(&asym_ids);
        if asym_ids.is_empty() {
            findings.push(malformed(category, row, "asym_id_list"));
            continue;
        }
        let Some(assembly) = assemblies.get_mut(assembly_id.as_ref()) else {
            findings.push(
                Diagnostic::new(Code::E6013)
                    .with_context("assembly", assembly_id)
                    .with_context("row", row.to_string()),
            );
            continue;
        };
        assembly.generators.push(Generator {
            oper_expression,
            asym_ids,
        });
    }
    if findings.is_empty() {
        Ok(())
    } else {
        Err(findings)
    }
}

fn operator_transform(
    category: &Category,
    row: usize,
    findings: &mut Vec<Diagnostic>,
) -> Option<Rigid> {
    let affine = match read_affine(category, row, Code::E6012) {
        Ok(transform) => transform,
        Err(finding) => {
            findings.push(finding);
            return None;
        }
    };
    let transform = Rigid::new(affine.matrix, affine.translation);
    if valid_rotation(&transform) {
        Some(transform)
    } else {
        findings.push(malformed(category, row, "matrix"));
        None
    }
}

fn valid_rotation(transform: &Rigid) -> bool {
    if !transform.translation.iter().all(|value| value.is_finite()) {
        return false;
    }
    let rows = transform.rotation;
    for left in 0..3 {
        for right in 0..3 {
            let dot = (0..3)
                .map(|axis| rows[left][axis] * rows[right][axis])
                .sum::<f64>();
            let expected = if left == right { 1.0 } else { 0.0 };
            if !dot.is_finite() || (dot - expected).abs() > 1e-4 {
                return false;
            }
        }
    }
    (transform.determinant().abs() - 1.0).abs() <= 1e-4
}

fn required_identifier(
    category: &Category,
    item: &str,
    row: usize,
    findings: &mut Vec<Diagnostic>,
) -> Option<Box<str>> {
    match category.identifier(item, row) {
        Some(value) if !value.trim().is_empty() => Some(value.into_owned().into_boxed_str()),
        _ => {
            findings.push(malformed(category, row, item));
            None
        }
    }
}

fn optional_text(category: &Category, item: &str, row: usize) -> Option<Box<str>> {
    category
        .identifier(item, row)
        .map(|value| value.into_owned().into_boxed_str())
}

fn optional_u32(
    category: &Category,
    item: &str,
    row: usize,
    findings: &mut Vec<Diagnostic>,
) -> Option<u32> {
    let value = category.value(item, row)?;
    if !value.is_recorded() {
        return None;
    }
    if let Some(value) = value
        .as_integer()
        .and_then(|value| u32::try_from(value).ok())
    {
        Some(value)
    } else {
        findings.push(malformed(category, row, item));
        None
    }
}

fn parse_asym_ids(text: &str) -> Box<[Box<str>]> {
    text.split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(Box::from)
        .collect()
}

fn malformed(category: &Category, row: usize, item: &str) -> Diagnostic {
    Diagnostic::new(Code::E6012)
        .with_context("category", category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

fn finish<T>(value: T, findings: Vec<Diagnostic>) -> Result<T, Vec<Diagnostic>> {
    if findings.is_empty() {
        Ok(value)
    } else {
        Err(findings)
    }
}

#[cfg(test)]
#[path = "lower_tests.rs"]
mod tests;
