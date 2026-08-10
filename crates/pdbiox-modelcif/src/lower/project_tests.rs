use super::lower;
use pdbiox_core::io::InputBuffer;

pub(crate) const MODEL_CIF: &str = "data_model\n\
loop_\n_ma_model_list.ordinal_id\n_ma_model_list.assembly_id\n_ma_model_list.model_name\n_ma_model_list.model_type\n\
1 1 prediction 'Ab initio model'\n#\n\
loop_\n_ma_target_entity.entity_id\n_ma_target_entity.data_id\n_ma_target_entity.origin\n1 1 'reference database'\n#\n\
loop_\n_ma_template_details.template_id\n_ma_template_details.target_asym_id\n_ma_template_details.template_name\n1 A template1\n#\n\
loop_\n_ma_protocol_step.protocol_id\n_ma_protocol_step.step_id\n_ma_protocol_step.method_type\n_ma_protocol_step.software_group_id\n1 1 modeling 1\n#\n\
loop_\n_ma_software_group.group_id\n_ma_software_group.software_id\n1 1\n#\n\
loop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n\
1 pLDDT pLDDT local\n2 PAE PAE local-pairwise\n3 pTM pTM global\n#\n\
loop_\n_ma_qa_metric_local.model_id\n_ma_qa_metric_local.label_asym_id\n_ma_qa_metric_local.label_seq_id\n_ma_qa_metric_local.label_comp_id\n_ma_qa_metric_local.metric_id\n_ma_qa_metric_local.metric_value\n\
1 A 1 ALA 1 91.5\n#\n\
loop_\n_ma_qa_metric_local_pairwise.model_id\n_ma_qa_metric_local_pairwise.label_asym_id_1\n_ma_qa_metric_local_pairwise.label_seq_id_1\n_ma_qa_metric_local_pairwise.label_asym_id_2\n_ma_qa_metric_local_pairwise.label_seq_id_2\n_ma_qa_metric_local_pairwise.metric_id\n_ma_qa_metric_local_pairwise.metric_value\n\
1 A 1 A 2 2 3.2\n#\n\
loop_\n_ma_qa_metric_global.model_id\n_ma_qa_metric_global.metric_id\n_ma_qa_metric_global.metric_value\n1 3 0.82\n#\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 0 0 0\n";

pub(crate) fn document() -> pdbiox_cif::Document {
    let input = InputBuffer::from_bytes(MODEL_CIF.as_bytes().to_vec());
    pdbiox_cif::parse(&input)
        .unwrap_or_else(|findings| panic!("fixture failed: {findings:?}"))
        .0
}

#[test]
fn lowers_metadata_and_three_confidence_modes() {
    let (model, findings) = lower(&document());
    assert!(findings.is_empty());
    assert_eq!(model.models.len(), 1);
    assert_eq!(model.targets.len(), 1);
    assert_eq!(model.templates.len(), 1);
    assert_eq!(model.protocol_steps.len(), 1);
    assert_eq!(model.software_groups.len(), 1);
    assert_eq!(model.confidence().plddt().count(), 1);
    assert_eq!(model.confidence().pae().count(), 1);
    assert_eq!(model.confidence().ptm().count(), 1);
    assert_eq!(model.categories.len(), 9);
}

#[test]
fn deprecated_model_identifier_is_reported_instead_of_substituted() {
    let source = "data_model\nloop_\n_ma_model_list.model_id\n_ma_model_list.model_name\n1 old\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (document, _) =
        pdbiox_cif::parse(&input).unwrap_or_else(|findings| panic!("fixture failed: {findings:?}"));
    let (model, findings) = lower(&document);
    assert!(model.models.is_empty());
    assert_eq!(findings.len(), 1);
    assert_eq!(findings[0].code(), pdbiox_core::Code::E2002);
}

#[test]
fn missing_data_block_is_reported_instead_of_treated_as_empty_metadata() {
    let document = pdbiox_cif::Document::new();
    let (model, findings) = lower(&document);
    assert!(model.categories.is_empty());
    assert_eq!(findings.len(), 1);
    assert_eq!(findings[0].code(), pdbiox_core::Code::E1106);
}

#[test]
fn ihm_categories_survive_the_lossless_document_round_trip() {
    let source = "data_ihm\nloop_\n_ihm_model_list.model_id\n_ihm_model_list.model_name\n1 'integrative model'\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (document, _) = pdbiox_cif::parse(&input)
        .unwrap_or_else(|findings| panic!("IHM fixture failed: {findings:?}"));
    let written = pdbiox_cif::write_preserving(&document);
    let reparsed = InputBuffer::from_bytes(written.into_bytes());
    let (round_trip, _) = pdbiox_cif::parse(&reparsed)
        .unwrap_or_else(|findings| panic!("IHM round trip failed: {findings:?}"));
    let category = round_trip
        .first_block()
        .and_then(|block| block.category("ihm_model_list"));
    assert_eq!(
        category.and_then(|value| value.text("model_name", 0)),
        Some("integrative model")
    );
}
