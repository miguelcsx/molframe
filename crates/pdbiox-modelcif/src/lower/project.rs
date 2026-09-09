use crate::{CompactColumn, MemoryBudget, ModelCategory, ModelCif, ModelCifError, ModelCifOptions};
use pdbiox_cif::Document;
use pdbiox_core::{Code, Diagnostic};

/// Explicitly materialises a CIF document's `ModelCIF` categories into one compact store.
///
/// Typed rows are validated without retaining a second representation. Unknown
/// categories and items remain available through the same columnar store.
///
/// # Errors
///
/// Returns a configuration, capacity, memory-ceiling, or allocation error
/// before returning a partial model.
pub fn lower(document: &Document) -> Result<(ModelCif, Vec<Diagnostic>), ModelCifError> {
    lower_with_options(document, ModelCifOptions::new())
}

/// Explicitly materialises `ModelCIF` categories under a memory policy.
///
/// # Errors
///
/// Returns a configuration, capacity, memory-ceiling, or allocation error
/// before returning a partial model.
pub fn lower_with_options(
    document: &Document,
    options: ModelCifOptions,
) -> Result<(ModelCif, Vec<Diagnostic>), ModelCifError> {
    let options = options.validate()?;
    let Some(block) = document.first_block() else {
        return Ok((ModelCif::default(), vec![Diagnostic::new(Code::E1106)]));
    };
    let count = block
        .categories()
        .filter(|category| category.name().starts_with("ma_"))
        .count();
    let mut categories = Vec::new();
    categories
        .try_reserve_exact(count)
        .map_err(|_| ModelCifError::Allocation)?;
    let mut budget = MemoryBudget::new(options.memory_limit_bytes);
    for category in block
        .categories()
        .filter(|category| category.name().starts_with("ma_"))
    {
        categories.push(ModelCategory::from_category(category, &mut budget)?);
    }
    let model = ModelCif { categories };
    let findings = validate_model(&model);
    Ok((model, findings))
}

pub(crate) fn validate_model(model: &ModelCif) -> Vec<Diagnostic> {
    let mut findings = Vec::new();
    validate(model.category("ma_model_list"), MODEL, &mut findings);
    validate(model.category("ma_target_entity"), TARGET, &mut findings);
    validate(
        model.category("ma_template_details"),
        TEMPLATE,
        &mut findings,
    );
    validate(model.category("ma_protocol_step"), PROTOCOL, &mut findings);
    validate(model.category("ma_software_group"), SOFTWARE, &mut findings);
    validate(model.category("ma_qa_metric"), DEFINITION, &mut findings);
    validate(model.category("ma_qa_metric_global"), GLOBAL, &mut findings);
    validate(model.category("ma_qa_metric_local"), LOCAL, &mut findings);
    validate(
        model.category("ma_qa_metric_local_pairwise"),
        PAIRWISE,
        &mut findings,
    );
    findings
}

#[derive(Clone, Copy)]
enum RequiredKind {
    Identifier,
    Integer,
    Number,
}

#[derive(Clone, Copy)]
struct Required {
    item: &'static str,
    kind: RequiredKind,
}

const fn id(item: &'static str) -> Required {
    Required {
        item,
        kind: RequiredKind::Identifier,
    }
}

const fn integer(item: &'static str) -> Required {
    Required {
        item,
        kind: RequiredKind::Integer,
    }
}

const fn number(item: &'static str) -> Required {
    Required {
        item,
        kind: RequiredKind::Number,
    }
}

const MODEL: &[Required] = &[id("ordinal_id")];
const TARGET: &[Required] = &[id("entity_id")];
const TEMPLATE: &[Required] = &[id("template_id")];
const PROTOCOL: &[Required] = &[id("protocol_id"), id("step_id"), id("method_type")];
const SOFTWARE: &[Required] = &[id("group_id"), id("software_id")];
const DEFINITION: &[Required] = &[id("id"), id("type"), id("mode")];
const GLOBAL: &[Required] = &[id("model_id"), id("metric_id"), number("metric_value")];
const LOCAL: &[Required] = &[
    id("model_id"),
    id("label_asym_id"),
    integer("label_seq_id"),
    id("metric_id"),
    number("metric_value"),
];
const PAIRWISE: &[Required] = &[
    id("model_id"),
    id("label_asym_id_1"),
    integer("label_seq_id_1"),
    id("label_asym_id_2"),
    integer("label_seq_id_2"),
    id("metric_id"),
    number("metric_value"),
];

fn validate(
    category: Option<&ModelCategory>,
    required: &[Required],
    findings: &mut Vec<Diagnostic>,
) {
    let Some(category) = category else {
        return;
    };
    let columns = required
        .iter()
        .map(|value| category.column(value.item))
        .collect::<Vec<_>>();
    for row in 0..category.row_count() {
        for (required, column) in required.iter().zip(&columns) {
            if let Err(finding) = validate_value(category, *required, *column, row) {
                findings.push(finding);
                break;
            }
        }
    }
}

fn validate_value(
    category: &ModelCategory,
    required: Required,
    column: Option<&CompactColumn>,
    row: usize,
) -> Result<(), Diagnostic> {
    let Some(column) = column else {
        return Err(missing_item(category, required.item, row));
    };
    match required.kind {
        RequiredKind::Identifier => column
            .identifier(row)
            .map(|_| ())
            .ok_or_else(|| missing_item(category, required.item, row)),
        RequiredKind::Integer => column
            .value(row)
            .ok_or_else(|| missing_item(category, required.item, row))?
            .integer()
            .map(|_| ())
            .ok_or_else(|| wrong_type(category, required.item, row)),
        RequiredKind::Number => column
            .value(row)
            .ok_or_else(|| missing_item(category, required.item, row))?
            .number()
            .filter(|value| value.is_finite())
            .map(|_| ())
            .ok_or_else(|| wrong_type(category, required.item, row)),
    }
}

fn missing_item(category: &ModelCategory, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2002)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

fn wrong_type(category: &ModelCategory, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2004)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

#[cfg(test)]
#[path = "project_tests.rs"]
pub(crate) mod tests;
