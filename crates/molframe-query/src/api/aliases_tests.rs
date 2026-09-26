use super::*;
use molframe_core::diagnostic::Code;

fn query(source: &str) -> Query {
    match Query::compile(source) {
        Ok(query) => query,
        Err(findings) => panic!("{source:?} did not compile: {findings:?}"),
    }
}

fn aliases(definitions: &[(&str, &str)]) -> QueryAliases {
    let mut aliases = QueryAliases::new();
    for (name, source) in definitions {
        if let Err(finding) = aliases.define(*name, query(source)) {
            panic!("{name} could not be defined: {finding:?}");
        }
    }
    aliases
}

fn error_code(result: Result<Query, Vec<Diagnostic>>) -> Code {
    match result {
        Ok(query) => panic!("resolution succeeded as {:?}", query.source()),
        Err(findings) => findings[0].code(),
    }
}

#[test]
fn a_defined_name_resolves_to_its_definition() {
    let aliases = aliases(&[("ligand", "resname ATP")]);
    let Ok(closed) = aliases.resolve(&query("$ligand")) else {
        panic!("ligand is defined");
    };
    assert_eq!(closed.fingerprint(), query("resname ATP").fingerprint());
    assert!(closed.references().is_empty());
}

#[test]
fn nested_names_resolve_through_every_level() {
    let aliases = aliases(&[
        ("heme", "resname HEM"),
        ("site", "within 5 of $heme"),
        ("pocket", "byres ($site) and protein"),
    ]);
    let Ok(closed) = aliases.resolve(&query("$pocket and not $heme")) else {
        panic!("every name is defined");
    };
    let expected = query("(byres (within 5 of (resname HEM)) and protein) and not (resname HEM)");
    assert_eq!(closed.fingerprint(), expected.fingerprint());
    assert_eq!(query(closed.source()).fingerprint(), closed.fingerprint());
}

#[test]
fn the_long_group_spelling_and_the_dollar_spelling_are_the_same_reference() {
    assert_eq!(
        query("group heme").fingerprint(),
        query("$heme").fingerprint()
    );
}

#[test]
fn an_unknown_name_is_reported_rather_than_matching_nothing() {
    let aliases = aliases(&[("heme", "resname HEM")]);
    assert_eq!(error_code(aliases.resolve(&query("$hem"))), Code::E4005);
}

#[test]
fn a_name_that_refers_to_itself_is_a_cycle() {
    let aliases = aliases(&[("loop", "protein or $loop")]);
    let Err(findings) = aliases.resolve(&query("$loop")) else {
        panic!("a self-reference cannot resolve");
    };
    assert_eq!(findings[0].code(), Code::E4006);
    assert!(
        findings[0]
            .context()
            .iter()
            .any(|item| item.value() == "loop -> loop")
    );
}

#[test]
fn a_cycle_through_several_names_is_reported_with_its_path() {
    let aliases = aliases(&[("a", "$b"), ("b", "$c"), ("c", "water or $a")]);
    let Err(findings) = aliases.resolve(&query("protein and $a")) else {
        panic!("a cycle cannot resolve");
    };
    assert_eq!(findings[0].code(), Code::E4006);
    assert!(
        findings[0]
            .context()
            .iter()
            .any(|item| item.value() == "a -> b -> c -> a")
    );
}

#[test]
fn a_reference_chain_deeper_than_the_bound_is_refused() {
    let mut definitions = Vec::new();
    for level in 0..=QueryAliases::MAX_DEPTH {
        definitions.push((format!("n{level}"), format!("$n{}", level + 1)));
    }
    definitions.push((
        format!("n{}", QueryAliases::MAX_DEPTH + 1),
        "all".to_owned(),
    ));
    let mut aliases = QueryAliases::new();
    for (name, source) in &definitions {
        assert!(aliases.define(name.as_str(), query(source)).is_ok());
    }
    assert_eq!(error_code(aliases.resolve(&query("$n0"))), Code::E4007);
}

#[test]
fn resolution_is_deterministic_and_independent_of_definition_order() {
    let forward = aliases(&[("a", "resname HEM"), ("b", "within 4 of $a")]);
    let backward = aliases(&[("b", "within 4 of $a"), ("a", "resname HEM")]);
    let target = query("$b or $a");
    let (Ok(left), Ok(right)) = (forward.resolve(&target), backward.resolve(&target)) else {
        panic!("both resolve");
    };
    assert_eq!(left.source(), right.source());
    assert_eq!(left.fingerprint(), right.fingerprint());
}

#[test]
fn redefining_a_dependency_changes_what_a_dependent_resolves_to() {
    let mut aliases = aliases(&[("ligand", "resname ATP"), ("pocket", "within 5 of $ligand")]);
    let pocket = query("$pocket");
    let Ok(before) = aliases.resolve(&pocket) else {
        panic!("resolves");
    };
    assert!(aliases.define("ligand", query("resname ADP")).is_ok());
    let Ok(after) = aliases.resolve(&pocket) else {
        panic!("resolves");
    };
    assert_ne!(before.fingerprint(), after.fingerprint());
    assert_eq!(
        after.fingerprint(),
        query("within 5 of (resname ADP)").fingerprint()
    );
}

#[test]
fn a_query_without_references_is_returned_unchanged() {
    let original = query("chain A and name CA");
    let Ok(resolved) = QueryAliases::new().resolve(&original) else {
        panic!("nothing to resolve");
    };
    assert_eq!(resolved.source(), original.source());
}

#[test]
fn references_are_sorted_and_unique() {
    assert_eq!(
        query("$b or ($a and within 3 of $b)").references(),
        vec!["a", "b"]
    );
}

#[test]
fn invalid_names_cannot_be_defined() {
    let mut aliases = QueryAliases::new();
    for name in ["", "1abc", "has space", "dash-ed", "$x"] {
        assert!(aliases.define(name, query("all")).is_err(), "{name:?}");
    }
    assert!(aliases.define("_ok_9", query("all")).is_ok());
}

#[test]
fn the_typed_builder_and_text_resolve_to_the_same_plan() {
    let built = Query::residues_within(5.0, query("$heme")) & Query::from(crate::col::is_protein());
    let aliases = aliases(&[("heme", "resname HEM")]);
    let (Ok(from_builder), Ok(from_text)) = (
        aliases.resolve(&built),
        aliases.resolve(&query("(byres (within 5 of ($heme))) and (protein)")),
    ) else {
        panic!("both resolve");
    };
    assert_eq!(from_builder.fingerprint(), from_text.fingerprint());
}
