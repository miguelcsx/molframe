use crate::{Groups, Query, col};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{AtomSelection, SecondaryStructure as Ss, Structure};
use std::fmt::Write as _;

#[test]
fn secondary_macros_select_exact_states_and_helpers_share_their_ir() {
    let states = [
        Ss::Unknown,
        Ss::Coil,
        Ss::AlphaHelix,
        Ss::ThreeTenHelix,
        Ss::PiHelix,
        Ss::OtherHelix,
        Ss::PolyProline,
        Ss::Strand,
        Ss::BetaBridge,
        Ss::Turn,
        Ss::Bend,
    ];
    let mut text = String::from(
        "data_secondary
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
",
    );
    for index in 1..=states.len() {
        writeln!(text, "ATOM {index} C CA GLY A {index} {index} 0 0").expect("write into String");
    }
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (raw, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let mut data = raw.data().clone();
    data.secondary_structure = states.to_vec().into();
    let structure = Structure::new(data);
    let cases = [
        (col::helix(), vec![2, 3, 4, 5, 6]),
        (col::strand(), vec![7]),
        (col::sheet(), vec![7, 8]),
        (col::alpha_helix(), vec![2]),
        (col::helix_310(), vec![3]),
        (col::pi_helix(), vec![4]),
        (col::polyproline(), vec![6]),
        (col::bridge(), vec![8]),
        (col::turn(), vec![9]),
        (col::bend(), vec![10]),
        (col::coil(), vec![1]),
    ];
    for (builder, expected) in cases {
        let source = builder.source().to_owned();
        let built = Query::from_builder(builder);
        let parsed = Query::compile(&source).expect("macro parses");
        assert_eq!(built.fingerprint(), parsed.fingerprint());
        let selected = parsed
            .evaluate_in(
                &structure,
                &AtomSelection::All(structure.atom_count()),
                &AnalysisPolicy::default(),
                &Groups::new(),
                None,
            )
            .expect("evaluate")
            .selection;
        assert_eq!(selected.iter().collect::<Vec<_>>(), expected, "{source}");
        let ast = crate::parser::parse(&source).expect("parse").expr;
        let printed = crate::print::print(&ast);
        assert_eq!(crate::parser::parse(&printed).expect("reparse").expr, ast);
        let completion = crate::complete(&source, source.len(), &crate::QueryAliases::new(), None);
        assert!(completion.items.iter().any(|item| item.label == source));
    }
}
