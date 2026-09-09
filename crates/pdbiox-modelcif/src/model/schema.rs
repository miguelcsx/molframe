pub(crate) fn identifier_item(category: &str, item: &str) -> bool {
    match category {
        "ma_model_list" => matches!(
            item,
            "ordinal_id" | "assembly_id" | "model_name" | "model_type"
        ),
        "ma_target_entity" => matches!(item, "entity_id" | "data_id" | "origin"),
        "ma_template_details" => matches!(
            item,
            "template_id"
                | "target_asym_id"
                | "template_name"
                | "template_origin"
                | "template_entity_type"
        ),
        "ma_protocol_step" => matches!(
            item,
            "protocol_id"
                | "step_id"
                | "method_type"
                | "step_name"
                | "details"
                | "software_group_id"
        ),
        "ma_software_group" => {
            matches!(item, "group_id" | "software_id" | "parameter_group_id")
        }
        "ma_qa_metric" => matches!(item, "id" | "name" | "type" | "mode" | "software_group_id"),
        "ma_qa_metric_global" => matches!(item, "model_id" | "metric_id"),
        "ma_qa_metric_local" => matches!(
            item,
            "model_id" | "label_asym_id" | "label_comp_id" | "metric_id"
        ),
        "ma_qa_metric_local_pairwise" => matches!(
            item,
            "model_id"
                | "label_asym_id_1"
                | "label_comp_id_1"
                | "label_asym_id_2"
                | "label_comp_id_2"
                | "metric_id"
        ),
        _ => false,
    }
}
