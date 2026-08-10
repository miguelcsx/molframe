use crate::{
    GlobalMetric, LocalMetric, MetricDefinition, ModelCategory, ModelCif, ModelDescription,
    PairwiseMetric, ProtocolStep, SoftwareGroup, Target, Template,
};
use pdbiox_cif::{Category, Document};
use pdbiox_core::{Code, Diagnostic};

/// Interprets every supported `ModelCIF` category and reports every malformed
/// typed row without discarding the lossless category representation.
#[must_use]
pub fn lower(document: &Document) -> (ModelCif, Vec<Diagnostic>) {
    let Some(block) = document.first_block() else {
        return (ModelCif::default(), vec![Diagnostic::new(Code::E1106)]);
    };
    let mut model = ModelCif {
        categories: block
            .categories()
            .filter(|category| category.name().starts_with("ma_"))
            .map(ModelCategory::from_category)
            .collect(),
        ..ModelCif::default()
    };
    let mut findings = Vec::new();
    if let Some(category) = block.category("ma_model_list") {
        lower_rows(
            category,
            model_description,
            &mut model.models,
            &mut findings,
        );
    }
    if let Some(category) = block.category("ma_target_entity") {
        lower_rows(category, target, &mut model.targets, &mut findings);
    }
    if let Some(category) = block.category("ma_template_details") {
        lower_rows(category, template, &mut model.templates, &mut findings);
    }
    if let Some(category) = block.category("ma_protocol_step") {
        lower_rows(category, protocol, &mut model.protocol_steps, &mut findings);
    }
    if let Some(category) = block.category("ma_software_group") {
        lower_rows(
            category,
            software,
            &mut model.software_groups,
            &mut findings,
        );
    }
    if let Some(category) = block.category("ma_qa_metric") {
        lower_rows(
            category,
            metric_definition,
            &mut model.quality.definitions,
            &mut findings,
        );
    }
    if let Some(category) = block.category("ma_qa_metric_global") {
        lower_rows(category, global, &mut model.quality.global, &mut findings);
    }
    if let Some(category) = block.category("ma_qa_metric_local") {
        lower_rows(category, local, &mut model.quality.local, &mut findings);
    }
    if let Some(category) = block.category("ma_qa_metric_local_pairwise") {
        lower_rows(
            category,
            pairwise,
            &mut model.quality.pairwise,
            &mut findings,
        );
    }
    (model, findings)
}

fn lower_rows<T>(
    category: &Category,
    parse: fn(&Category, usize) -> Result<T, Diagnostic>,
    output: &mut Vec<T>,
    findings: &mut Vec<Diagnostic>,
) {
    for row in 0..category.row_count() {
        match parse(category, row) {
            Ok(value) => output.push(value),
            Err(finding) => findings.push(finding),
        }
    }
}

fn id(category: &Category, item: &str, row: usize) -> Option<Box<str>> {
    category
        .identifier(item, row)
        .map(|value| value.into_owned().into_boxed_str())
}

fn model_description(category: &Category, row: usize) -> Result<ModelDescription, Diagnostic> {
    Ok(ModelDescription {
        id: required_id(category, "ordinal_id", row)?,
        assembly_id: id(category, "assembly_id", row),
        name: id(category, "model_name", row),
        model_type: id(category, "model_type", row),
    })
}

fn target(category: &Category, row: usize) -> Result<Target, Diagnostic> {
    Ok(Target {
        entity_id: required_id(category, "entity_id", row)?,
        data_id: id(category, "data_id", row),
        origin: id(category, "origin", row),
    })
}

fn template(category: &Category, row: usize) -> Result<Template, Diagnostic> {
    Ok(Template {
        id: required_id(category, "template_id", row)?,
        target_chain_id: id(category, "target_asym_id", row),
        name: id(category, "template_name", row),
        origin: id(category, "template_origin", row),
        entity_type: id(category, "template_entity_type", row),
    })
}

fn protocol(category: &Category, row: usize) -> Result<ProtocolStep, Diagnostic> {
    Ok(ProtocolStep {
        protocol_id: required_id(category, "protocol_id", row)?,
        step_id: required_id(category, "step_id", row)?,
        method_type: required_id(category, "method_type", row)?,
        name: id(category, "step_name", row),
        details: id(category, "details", row),
        software_group_id: id(category, "software_group_id", row),
    })
}

fn software(category: &Category, row: usize) -> Result<SoftwareGroup, Diagnostic> {
    Ok(SoftwareGroup {
        group_id: required_id(category, "group_id", row)?,
        software_id: required_id(category, "software_id", row)?,
        parameter_group_id: id(category, "parameter_group_id", row),
    })
}

fn metric_definition(category: &Category, row: usize) -> Result<MetricDefinition, Diagnostic> {
    Ok(MetricDefinition {
        id: required_id(category, "id", row)?,
        name: id(category, "name", row),
        metric_type: required_id(category, "type", row)?,
        mode: required_id(category, "mode", row)?,
        software_group_id: id(category, "software_group_id", row),
    })
}

fn global(category: &Category, row: usize) -> Result<GlobalMetric, Diagnostic> {
    Ok(GlobalMetric {
        model_id: required_id(category, "model_id", row)?,
        metric_id: required_id(category, "metric_id", row)?,
        value: required_number(category, "metric_value", row)?,
    })
}

fn local(category: &Category, row: usize) -> Result<LocalMetric, Diagnostic> {
    Ok(LocalMetric {
        model_id: required_id(category, "model_id", row)?,
        chain_id: required_id(category, "label_asym_id", row)?,
        sequence_id: required_integer(category, "label_seq_id", row)?,
        component_id: id(category, "label_comp_id", row),
        metric_id: required_id(category, "metric_id", row)?,
        value: required_number(category, "metric_value", row)?,
    })
}

fn pairwise(category: &Category, row: usize) -> Result<PairwiseMetric, Diagnostic> {
    Ok(PairwiseMetric {
        model_id: required_id(category, "model_id", row)?,
        first_chain_id: required_id(category, "label_asym_id_1", row)?,
        first_sequence_id: required_integer(category, "label_seq_id_1", row)?,
        second_chain_id: required_id(category, "label_asym_id_2", row)?,
        second_sequence_id: required_integer(category, "label_seq_id_2", row)?,
        metric_id: required_id(category, "metric_id", row)?,
        value: required_number(category, "metric_value", row)?,
    })
}

fn required_id(category: &Category, item: &str, row: usize) -> Result<Box<str>, Diagnostic> {
    id(category, item, row).ok_or_else(|| missing_item(category, item, row))
}

fn required_number(category: &Category, item: &str, row: usize) -> Result<f64, Diagnostic> {
    let value = category
        .value(item, row)
        .ok_or_else(|| missing_item(category, item, row))?;
    value
        .as_float()
        .filter(|number| number.is_finite())
        .ok_or_else(|| wrong_type(category, item, row))
}

fn required_integer(category: &Category, item: &str, row: usize) -> Result<i64, Diagnostic> {
    let value = category
        .value(item, row)
        .ok_or_else(|| missing_item(category, item, row))?;
    value
        .as_integer()
        .ok_or_else(|| wrong_type(category, item, row))
}

fn missing_item(category: &Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2002)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

fn wrong_type(category: &Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2004)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

#[cfg(test)]
#[path = "project_tests.rs"]
pub(crate) mod tests;
