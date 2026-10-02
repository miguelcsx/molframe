use super::*;

fn parsed(source: &str) -> Expr {
    match crate::parser::parse(source) {
        Ok(parsed) => parsed.expr,
        Err(findings) => panic!("{source:?} did not parse: {findings:?}"),
    }
}

const CORPUS: &[&str] = &[
    "all",
    "none",
    "protein and not water",
    "resname HEM or (chain A and name CA CB)",
    "byres (within 5 of resname HEM) and protein",
    "same residue as (name FE)",
    "same resname as (index 3)",
    "bonded 2 (element Fe)",
    "beyond 3.5 of ligand",
    "around 4 hetero",
    "sphzone 6 (resname ATP)",
    "sphlayer 2 6.5 (resname ATP)",
    "isolayer 1 3 (resname ATP)",
    "cyzone 5 2 -2 (protein)",
    "cylayer 1 5 2 -2 (protein)",
    "point 1.5 -2 3 4.25",
    "bfactor > 30.5 and occupancy <= 0.5",
    "prop abs charge >= 0.25",
    "10 < bfactor",
    "name 'C A' \"O'1\" and resname \\and",
    "name C* and not name CA",
    "atom A 12 CA",
    "global (chain B)",
    "$heme or group \"odd name\"",
    "chirality R",
    "backbone and nucleicbase or nucleicsugar",
    "altloc A and icode B and segid S1",
    "assembly 0:3",
    "same assembly as (index 0)",
];

#[test]
fn printed_queries_parse_back_to_the_same_plan() {
    for source in CORPUS {
        let expr = parsed(source);
        let printed = print(&expr);
        let reparsed = parsed(&printed);
        assert_eq!(expr, reparsed, "{source:?} printed as {printed:?}");
    }
}

#[test]
fn printing_is_idempotent() {
    for source in CORPUS {
        let once = print(&parsed(source));
        let twice = print(&parsed(&once));
        assert_eq!(once, twice, "{source:?}");
    }
}

#[test]
fn membership_values_that_the_lexer_would_split_survive_printing() {
    for value in ["C A", "a\"b", "back\\slash", "\\leading", "(x)", "OR", ""] {
        let expr = Expr::Membership {
            column: Column::AtomName,
            values: vec![value.into()],
        };
        let printed = print(&expr);
        if value.is_empty() {
            // An empty membership pattern still prints as a quoted value.
            assert!(printed.ends_with("\"\""), "{printed}");
            continue;
        }
        assert_eq!(parsed(&printed), expr, "{value:?} printed as {printed:?}");
    }
}

#[test]
fn a_parenthesised_print_never_triggers_the_mixed_precedence_warning() {
    let printed = print(&parsed("all or none and all"));
    let Ok(reparsed) = crate::parser::parse(&printed) else {
        panic!("printed query parses");
    };
    assert!(reparsed.warnings.is_empty(), "{:?}", reparsed.warnings);
}
