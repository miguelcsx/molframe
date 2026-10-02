use super::*;
use molframe_bench::{Sample, structure, structure_with_atoms};
use molframe_core::contract::AnalysisPolicy;

/// The atoms whose element is one of `symbols`, found by reading every atom.
fn brute_force(structure: &Structure, symbols: &[&str]) -> Vec<u32> {
    let wanted: Vec<Element> = symbols
        .iter()
        .filter_map(|symbol| Element::from_symbol(symbol))
        .collect();
    structure
        .data()
        .atoms()
        .filter(|atom| {
            atom.element()
                .is_some_and(|element| wanted.contains(&element))
        })
        .map(|atom| atom.index().get())
        .collect()
}

fn evaluated(structure: &Structure, source: &str) -> Vec<u32> {
    let query = crate::Query::compile(source).expect("the corpus compiles");
    let evaluation = query
        .evaluate(
            structure,
            &AnalysisPolicy::default(),
            &crate::Groups::new(),
            None,
        )
        .expect("the corpus evaluates");
    evaluation.selection.iter().collect()
}

#[test]
fn pruned_element_queries_equal_a_full_read_of_every_atom() {
    let corpus: [&[&str]; 10] = [
        &["C"],
        &["N"],
        &["O"],
        &["S"],
        &["FE"],
        &["Fe"],
        &["ZN"],
        &["C", "N"],
        &["FE", "ZN"],
        &["H", "MG", "ZN"],
    ];
    for source in [structure(Sample::Medium), structure(Sample::Large)] {
        for symbols in corpus {
            let text = format!("element {}", symbols.join(" "));
            assert_eq!(
                evaluated(&source, &text),
                brute_force(&source, symbols),
                "{text}"
            );
        }
    }
    let tiled = structure_with_atoms(100_000);
    for symbols in [&["C"][..], &["ZN"], &["S", "ZN"]] {
        let text = format!("element {}", symbols.join(" "));
        assert_eq!(
            evaluated(&tiled, &text),
            brute_force(&tiled, symbols),
            "{text}"
        );
    }
}

#[test]
fn an_element_the_summaries_exclude_admits_no_chunk() {
    let tiled = structure_with_atoms(100_000);
    assert!(tiled.data().chunks.len() > 1, "the fixture spans chunks");
    assert!(chunks_with_any_element(&tiled, &[Element::ZINC]).is_empty());
    let carbon = chunks_with_any_element(&tiled, &[Element::CARBON]);
    assert_eq!(carbon.len(), u64::from(tiled.atom_count()));
}

#[test]
fn a_chunk_holding_the_element_is_admitted_whole() {
    let medium = structure(Sample::Medium);
    let admitted = chunks_with_any_element(&medium, &[Element::IRON]);
    let iron = brute_force(&medium, &["FE"]);
    assert!(!iron.is_empty(), "4HHB carries haem iron");
    assert!(iron.iter().all(|atom| admitted.contains(*atom)));
}
