//! The heme pocket of 4HHB, written with named queries and without.
//!
//! The same pocket is selected twice: once as `$pocket` resolved through
//! definitions that refer to `$heme`, and once as the equivalent closed text.
//! The two must select identical atoms, and redefining `heme` must move the
//! pocket with it.

use molframe::{AnalysisPolicy, Query, QueryAliases, QueryStructure};
use molframe_bench::{Sample, structure};

fn compile(source: &str) -> Query {
    match Query::compile(source) {
        Ok(query) => query,
        Err(findings) => panic!("{source:?}: {findings:?}"),
    }
}

fn atoms(structure: &molframe_core::Structure, query: &Query) -> Vec<u32> {
    match structure.select_query(query, &AnalysisPolicy::default()) {
        Ok(evaluation) => evaluation.selection.iter().collect(),
        Err(findings) => panic!("{:?}: {findings:?}", query.source()),
    }
}

#[test]
fn the_named_heme_pocket_selects_the_same_atoms_as_its_closed_text() {
    let structure = structure(Sample::Medium);
    let mut aliases = QueryAliases::new();
    assert!(aliases.define("heme", compile("resname HEM")).is_ok());
    assert!(
        aliases
            .define("pocket", compile("byres (within 5 of $heme) and protein"))
            .is_ok()
    );
    let Ok(pocket) = aliases.resolve(&compile("$pocket")) else {
        panic!("every name is defined");
    };
    let closed = compile("byres (within 5 of (resname HEM)) and protein");
    let named = atoms(&structure, &pocket);
    assert_eq!(named, atoms(&structure, &closed));
    assert_eq!(pocket.fingerprint(), closed.fingerprint());
    // Four hemes of 43 atoms each.
    assert_eq!(atoms(&structure, &compile("resname HEM")).len(), 172);
    assert!(!named.is_empty());

    assert!(
        aliases
            .define("heme", compile("resname HEM and chain A"))
            .is_ok()
    );
    let Ok(moved) = aliases.resolve(&compile("$pocket")) else {
        panic!("every name is defined");
    };
    let narrowed = atoms(&structure, &moved);
    assert!(!narrowed.is_empty());
    assert!(narrowed.len() < named.len());
    assert!(
        narrowed
            .iter()
            .all(|atom| named.binary_search(atom).is_ok())
    );
}
