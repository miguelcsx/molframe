//! Enforces the promotion rules of `site/content/docs/parity.mdx`.
//!
//! A row at `✓` must name a golden workflow, and every golden workflow the
//! ledger names must have a `gw_NNN_` test.

const LEDGER: &str = include_str!("../../../site/content/docs/parity.mdx");

/// Every test source that may hold a golden workflow.
const GOLDEN_SOURCES: &[&str] = &[
    include_str!("golden_tests.rs"),
    include_str!("golden_workflows.rs"),
    include_str!("golden/connection_round_trip.rs"),
    include_str!("golden/named_queries.rs"),
    include_str!("golden/native_pending.rs"),
    include_str!("golden/native_pending/analysis.rs"),
    include_str!("golden/native_pending/compare_interop.rs"),
    include_str!("golden/native_pending/trajectory.rs"),
    include_str!("golden/native_pending/xtal.rs"),
];

const PROMOTED: &str = "✓";

struct Row<'a> {
    capability: &'a str,
    status: &'a str,
    workflow: &'a str,
}

/// Table body rows: `| capability | status | crate | GW | evidence |`.
fn rows(text: &str) -> Vec<Row<'_>> {
    text.lines()
        .filter(|line| line.starts_with('|'))
        .filter_map(|line| {
            let cells: Vec<&str> = line.split('|').map(str::trim).collect();
            let [_, capability, status, _, workflow, ..] = cells.as_slice() else {
                return None;
            };
            let is_header = *capability == "Capability" || capability.starts_with("---");
            (!is_header).then_some(Row {
                capability,
                status,
                workflow,
            })
        })
        .collect()
}

/// Three-digit workflow numbers in a cell, written `GW-NNN`.
fn workflow_numbers(cell: &str) -> Vec<&str> {
    cell.match_indices("GW-")
        .filter_map(|(start, _)| cell.get(start + 3..start + 6))
        .filter(|digits| digits.bytes().all(|byte| byte.is_ascii_digit()))
        .collect()
}

fn promoted_rows_without_a_workflow(text: &str) -> Vec<String> {
    rows(text)
        .into_iter()
        .filter(|row| row.status == PROMOTED && workflow_numbers(row.workflow).is_empty())
        .map(|row| row.capability.to_owned())
        .collect()
}

fn named_workflows_without_a_test(text: &str, sources: &[&str]) -> Vec<String> {
    let mut missing = Vec::new();
    for row in rows(text) {
        for number in workflow_numbers(row.workflow) {
            let needle = format!("fn gw_{number}_");
            if !sources.iter().any(|source| source.contains(&needle))
                && !missing.contains(&number.to_owned())
            {
                missing.push(number.to_owned());
            }
        }
    }
    missing
}

#[test]
fn the_ledger_lists_the_whole_register() {
    assert!(
        rows(LEDGER).len() >= 300,
        "only {} rows parsed; the table layout changed",
        rows(LEDGER).len()
    );
}

#[test]
fn every_promoted_row_names_a_golden_workflow() {
    let bare = promoted_rows_without_a_workflow(LEDGER);
    assert!(bare.is_empty(), "promoted without a GW: {bare:?}");
}

#[test]
fn every_named_workflow_has_a_test() {
    let missing = named_workflows_without_a_test(LEDGER, GOLDEN_SOURCES);
    assert!(
        missing.is_empty(),
        "no `fn gw_NNN_` test for GW-{missing:?}"
    );
}

#[test]
fn the_checks_reject_a_bare_promotion_and_a_ghost_workflow() {
    let ledger = "| Capability | molframe | Crate | GW | Evidence |\n\
                  | --- | :--: | --- | --- | --- |\n\
                  | Bare | ✓ | `molframe-core` | — | none |\n\
                  | Named | ✓ | `molframe-core` | GW-999 | gemmi 0.7.5 |\n";
    assert_eq!(promoted_rows_without_a_workflow(ledger), ["Bare"]);
    assert_eq!(named_workflows_without_a_test(ledger, &[]), ["999"]);
    assert!(named_workflows_without_a_test(ledger, &["fn gw_999_x() {}"]).is_empty());
}
