use super::{NcsBindOptions, bind_ncs, expand_ncs};
use crate::view::tests::attached;
use crate::{NcsDomains, NcsSet, lower_ncs};
use molframe_cif::parse;
use molframe_core::{Code, InputBuffer};

const ATOMS: &str = r"data_ncs
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 1 A 0 0 0
2 C CA ALA A 2 2 A 3 0 0
3 C CA SER A 3 3 A 3 3 1
4 C CA GLY B 1 1 B 10 0 0
5 C CA ALA B 2 2 B 13 0 0
6 C CA SER B 3 3 B 13 3 1
7 C CA TRP C 1 1 C 50 50 50
8 C CA TRP C 2 2 C 53 50 50
9 C CA TRP C 3 3 C 53 53 51
";

const OPERATORS: &str = r"loop_
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
2 generate 1 0 0 0 0 1 0 0 0 0 1 100
";

const DECLARED: &str = r"loop_
_struct_ncs_ens.id
ens_1
loop_
_struct_ncs_dom.id
_struct_ncs_dom.pdbx_ens_id
d1 ens_1
d2 ens_1
loop_
_struct_ncs_dom_lim.dom_id
_struct_ncs_dom_lim.pdbx_ens_id
_struct_ncs_dom_lim.beg_label_asym_id
_struct_ncs_dom_lim.beg_label_seq_id
_struct_ncs_dom_lim.end_label_asym_id
_struct_ncs_dom_lim.end_label_seq_id
d1 ens_1 A 1 A 3
d2 ens_1 B 1 B 3
loop_
_struct_ncs_ens_gen.ens_id
_struct_ncs_ens_gen.dom_id_1
_struct_ncs_ens_gen.dom_id_2
_struct_ncs_ens_gen.oper_id
ens_1 d1 d2 1
";

fn set(extra: &str) -> NcsSet {
    let input = InputBuffer::from_bytes(format!("{ATOMS}{OPERATORS}{extra}").into_bytes());
    match parse(&input) {
        Ok((document, _)) => lower_ncs(&document).unwrap_or_else(|f| panic!("lowering: {f:?}")),
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    }
}

fn labels(structure: &molframe_core::Structure, pairs: &[super::NcsChainPair]) -> Vec<String> {
    let label = |chain| {
        structure
            .data()
            .chain(chain)
            .and_then(molframe_core::structure::ChainRef::label)
            .unwrap_or("?")
            .to_owned()
    };
    pairs
        .iter()
        .map(|pair| format!("{}>{}", label(pair.source), label(pair.target)))
        .collect()
}

#[test]
fn geometry_binds_the_given_operator_to_the_chain_it_maps() {
    let structure = attached(ATOMS);
    let pairs = bind_ncs(&structure, &set(""), &NcsBindOptions::default())
        .unwrap_or_else(|f| panic!("bind failed: {f}"));
    assert_eq!(labels(&structure, &pairs), ["A>B"]);
    assert!(pairs[0].rmsd.is_some_and(|value| value < 1e-4));
    assert!(!pairs[0].declared);
    assert_eq!(pairs[0].atoms, 3);
}

#[test]
fn a_tight_tolerance_rejects_a_pair_that_does_not_superpose() {
    let structure = attached(&ATOMS.replace("13 3 1", "13 3 9"));
    let pairs = bind_ncs(&structure, &set(""), &NcsBindOptions::default())
        .unwrap_or_else(|f| panic!("bind failed: {f}"));
    assert!(pairs.is_empty());
}

#[test]
fn declared_pairings_bind_by_domain_range_and_report_the_fit() {
    let structure = attached(ATOMS);
    let pairs = bind_ncs(&structure, &set(DECLARED), &NcsBindOptions::default())
        .unwrap_or_else(|f| panic!("bind failed: {f}"));
    assert_eq!(labels(&structure, &pairs), ["A>B"]);
    assert!(pairs[0].declared);
    assert!(pairs[0].rmsd.is_some_and(|value| value < 1e-4));
}

#[test]
fn a_domain_naming_an_absent_chain_is_reported() {
    let structure = attached(ATOMS);
    let bad = DECLARED.replace("B 1 B 3", "Z 1 Z 3");
    let Err(finding) = bind_ncs(&structure, &set(&bad), &NcsBindOptions::default()) else {
        panic!("absent chain accepted")
    };
    assert_eq!(finding.code(), Code::E6014);
}

#[test]
fn the_search_is_bounded_by_the_pair_limit() {
    let structure = attached(ATOMS);
    let options = NcsBindOptions {
        pair_limit: 1,
        ..NcsBindOptions::default()
    };
    let Err(finding) = bind_ncs(&structure, &set(""), &options) else {
        panic!("limit ignored")
    };
    assert_eq!(finding.code(), Code::E6015);
}

#[test]
fn expansion_appends_one_transformed_copy_per_generate_operator() {
    let structure = attached(ATOMS);
    let expanded =
        expand_ncs(&structure, &set("")).unwrap_or_else(|f| panic!("expand failed: {f:?}"));
    assert_eq!(expanded.atom_count(), structure.atom_count() * 2);
    let source = structure.positions();
    let moved = &expanded.positions()[source.len()..];
    for (before, after) in source.iter().zip(moved) {
        assert!((after[2] - before[2] - 100.0).abs() < 1e-4);
        assert!((after[0] - before[0]).abs() < 1e-4);
    }
    assert_eq!(&expanded.positions()[..source.len()], source);
}

#[test]
fn expansion_without_generate_operators_returns_the_structure() {
    let structure = attached(ATOMS);
    let none = NcsSet::from_operators([]).and_then(|set| set.with_domains(NcsDomains::default()));
    let Ok(none) = none else {
        panic!("empty set rejected")
    };
    let expanded = expand_ncs(&structure, &none).unwrap_or_else(|f| panic!("{f:?}"));
    assert_eq!(expanded.atom_count(), structure.atom_count());
}
