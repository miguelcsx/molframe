use super::{NCS_EXTENSION, NcsCode, NcsExt, lower_ncs};
use crate::view::tests::{ENTRY, attached};
use pdbiox_cif::parse;
use pdbiox_core::{Code, InputBuffer, ModelIndex};

const NCS: &str = r"data_ncs
loop_
_struct_ncs_oper.id
_struct_ncs_oper.code
_struct_ncs_oper.details
_struct_ncs_oper.matrix[1][1]
_struct_ncs_oper.matrix[1][2]
_struct_ncs_oper.matrix[1][3]
_struct_ncs_oper.vector[1]
_struct_ncs_oper.matrix[2][1]
_struct_ncs_oper.matrix[2][2]
_struct_ncs_oper.matrix[2][3]
_struct_ncs_oper.vector[2]
_struct_ncs_oper.matrix[3][1]
_struct_ncs_oper.matrix[3][2]
_struct_ncs_oper.matrix[3][3]
_struct_ncs_oper.vector[3]
given-op given 'already deposited' 1 0.2 0 0 0 1 0 0 0 0 1 0
new-op generate 'missing copy' 1 0 0 5 0 1 0 6 0 0 1 7
";

#[test]
fn lowering_preserves_non_rigid_given_operators_and_generate_semantics() {
    let set = parsed_set(NCS);
    assert_eq!(set.len(), 2);
    assert_eq!(set.generators().count(), 1);
    assert_eq!(
        set.get("given-op").map(|operator| operator.code),
        Some(NcsCode::Given)
    );
    assert_eq!(
        set.get("given-op")
            .map(|operator| operator.transform.matrix[0][1]),
        Some(0.2)
    );
}

#[test]
fn generated_view_applies_only_generate_operators_without_copying_source_storage() {
    let structure = attached(ENTRY).with_extension(NCS_EXTENSION, parsed_set(NCS));
    let Some(view) = structure.ncs_generated() else {
        panic!("NCS extension absent")
    };
    assert_eq!(view.copy_count(), 1);
    assert_eq!(view.atoms().count(), structure.atom_count() as usize);
    let first = view
        .positions(ModelIndex::new(0))
        .next()
        .and_then(|(_, position)| position);
    assert_eq!(first, Some([6.0, 6.0, 7.0]));
    assert!(
        structure.positions()[0]
            .iter()
            .zip([1.0, 0.0, 0.0])
            .all(|(actual, expected)| (*actual - expected).abs() < 1e-6)
    );
}

#[test]
fn an_unknown_controlled_code_is_rejected() {
    let invalid = NCS.replace("given-op given", "given-op maybe");
    let input = InputBuffer::from_bytes(invalid.into_bytes());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    };
    let Err(findings) = lower_ncs(&document) else {
        panic!("unknown NCS code accepted")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E6014));
}

fn parsed_set(text: &str) -> super::NcsSet {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    };
    match lower_ncs(&document) {
        Ok(set) => set,
        Err(findings) => panic!("NCS lowering failed: {findings:?}"),
    }
}
