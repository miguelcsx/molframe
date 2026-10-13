use crate::AffineTransform;
use crate::{NcsDomains, NcsOperator, NcsSet, lower_ncs, ncs_cif};
use molframe_cif::parse;
use molframe_core::{Code, InputBuffer};

const TEXT: &str = r"data_dom
loop_
_struct_ncs_oper.id
_struct_ncs_oper.code
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
1 given 1 0 0 10 0 1 0 0 0 0 1 0
loop_
_struct_ncs_ens.id
_struct_ncs_ens.point_group
ens_1 C2
loop_
_struct_ncs_dom.id
_struct_ncs_dom.pdbx_ens_id
d1 ens_1
d2 ens_1
loop_
_struct_ncs_dom_lim.dom_id
_struct_ncs_dom_lim.pdbx_ens_id
_struct_ncs_dom_lim.pdbx_component_id
_struct_ncs_dom_lim.beg_label_asym_id
_struct_ncs_dom_lim.beg_label_comp_id
_struct_ncs_dom_lim.beg_label_seq_id
_struct_ncs_dom_lim.end_label_asym_id
_struct_ncs_dom_lim.end_label_comp_id
_struct_ncs_dom_lim.end_label_seq_id
d1 ens_1 1 A GLY 1 A SER 3
d1 ens_1 2 A GLY 7 A SER 9
d2 ens_1 1 B GLY 1 B SER 3
loop_
_struct_ncs_ens_gen.ens_id
_struct_ncs_ens_gen.dom_id_1
_struct_ncs_ens_gen.dom_id_2
_struct_ncs_ens_gen.oper_id
ens_1 d1 d2 1
";

fn lowered(text: &str) -> Result<NcsSet, Vec<molframe_core::Diagnostic>> {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match parse(&input) {
        Ok((document, _)) => lower_ncs(&document),
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    }
}

#[test]
fn ensembles_domains_limits_and_pairings_are_lowered_in_order() {
    let set = lowered(TEXT).unwrap_or_else(|f| panic!("lowering: {f:?}"));
    let domains = set.domains();
    assert_eq!(domains.ensembles.len(), 1);
    assert_eq!(domains.ensembles[0].point_group.as_deref(), Some("C2"));
    assert_eq!(domains.domains.len(), 2);
    let first = domains.domain("d1").unwrap_or_else(|| panic!("d1 absent"));
    assert_eq!(first.limits.len(), 2);
    assert_eq!(first.limits[1].begin.seq_id, Some(7));
    assert_eq!(&*first.limits[0].end.asym_id, "A");
    assert_eq!(domains.pairings.len(), 1);
    assert_eq!(&*domains.pairings[0].operator, "1");
}

#[test]
fn written_domains_read_back_identically() {
    let set = lowered(TEXT).unwrap_or_else(|f| panic!("lowering: {f:?}"));
    let text = format!("data_again\n{}", ncs_cif(&set));
    let again = lowered(&text).unwrap_or_else(|f| panic!("re-read: {f:?}"));
    assert_eq!(again.domains(), set.domains());
    assert_eq!(again.len(), set.len());
}

#[test]
fn a_pairing_naming_an_unknown_operator_is_rejected_on_read() {
    let findings = lowered(&TEXT.replace("ens_1 d1 d2 1", "ens_1 d1 d2 9"))
        .err()
        .unwrap_or_else(|| panic!("unknown operator accepted"));
    assert!(findings.iter().any(|finding| finding.code() == Code::E6014));
}

#[test]
fn a_pairing_naming_an_unknown_domain_is_rejected_when_attached() {
    let operator = NcsOperator {
        id: "1".into(),
        code: crate::NcsCode::Given,
        details: None,
        transform: AffineTransform::IDENTITY,
    };
    let set = NcsSet::from_operators([operator]).unwrap_or_else(|f| panic!("{f}"));
    let mut domains = NcsDomains::default();
    domains.pairings.push(crate::NcsPairing {
        ensemble: "e".into(),
        domain_1: "x".into(),
        domain_2: "y".into(),
        operator: "1".into(),
    });
    let rejected = set.with_domains(domains).err();
    assert_eq!(rejected.map(|finding| finding.code()), Some(Code::E6014));
}

#[test]
fn a_limit_for_an_undeclared_domain_is_rejected() {
    let findings = lowered(&TEXT.replace("d2 ens_1 1 B", "d9 ens_1 1 B"))
        .err()
        .unwrap_or_else(|| panic!("orphan limit accepted"));
    assert!(findings.iter().any(|finding| finding.code() == Code::E6014));
}
