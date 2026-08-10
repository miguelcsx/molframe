use crate::lower;
use pdbiox_core::io::{InputBuffer, ReadOptions};

#[test]
fn canonical_write_keeps_modelcif_categories_and_metrics() {
    let document = crate::lower::tests::document();
    let (structure, _) = pdbiox_cif::lower(&document, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("lower failed: {findings:?}"));
    let (model, findings) = lower(&document);
    assert!(findings.is_empty());
    let options = pdbiox_cif::CifWriteOptions::new()
        .with_block_id("model")
        .with_generated_connection_ids()
        .with_connection_type_id("covale");
    let text = super::write_canonical_with_options(&structure, &model, &options)
        .unwrap_or_else(|error| panic!("write failed: {error}"));
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (round_trip, _) = pdbiox_cif::parse(&input)
        .unwrap_or_else(|findings| panic!("written `ModelCIF` failed: {findings:?}"));
    let (model, findings) = lower(&round_trip);
    assert!(findings.is_empty());
    assert_eq!(model.confidence().plddt().count(), 1);
    assert_eq!(model.confidence().pae().count(), 1);
    assert_eq!(model.models[0].name.as_deref(), Some("prediction"));
}
