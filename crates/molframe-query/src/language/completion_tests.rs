use super::*;
use crate::Query;

#[test]
fn completion_uses_byte_spans_for_unicode_prefixes() {
    let aliases = QueryAliases::new();
    let source = "α pr";
    let result = complete(source, source.len(), &aliases, None);
    assert_eq!(
        &source[result.replacement_start..result.replacement_end],
        "pr"
    );
    assert!(result.items.iter().any(|item| item.label == "protein"));
}

#[test]
fn completion_exposes_only_defined_aliases_after_dollar() {
    let mut aliases = QueryAliases::new();
    let query = Query::compile("protein").expect("query compiles");
    aliases
        .define("pocket", query)
        .expect("alias name is valid");
    let result = complete("$poc", 4, &aliases, None);
    assert_eq!(result.items.len(), 1);
    assert_eq!(result.items[0].label, "pocket");
    assert_eq!(result.items[0].kind, CompletionKind::Alias);
}

#[test]
fn completion_rejects_cursor_inside_utf8_scalar() {
    let aliases = QueryAliases::new();
    let result = complete("α", 1, &aliases, None);
    assert!(result.items.is_empty());
    assert_eq!(result.diagnostics.len(), 1);
}

struct Chains<'a> {
    labels: &'a [&'a str],
}
impl StructureValues for Chains<'_> {
    fn chain_labels(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        Box::new(self.labels.iter().copied())
    }
}

fn complete_with_values(
    source: &str,
    cursor: usize,
    aliases: &QueryAliases,
    values: Option<&dyn StructureValues>,
) -> CompletionResult {
    complete(source, cursor, aliases, values)
}

#[test]
fn completion_offers_structure_values_only_after_value_columns() {
    let aliases = QueryAliases::new();
    let chains = Chains {
        labels: &["A", "B"],
    };
    let after_chain = complete_with_values("chain ", 6, &aliases, Some(&chains));
    assert_eq!(
        after_chain
            .items
            .iter()
            .map(|item| (item.label.as_str(), item.kind))
            .collect::<Vec<_>>(),
        [("A", CompletionKind::Value), ("B", CompletionKind::Value)]
    );

    let mid = complete_with_values("chain A", 7, &aliases, Some(&chains));
    assert_eq!(&"chain A"[mid.replacement_start..mid.replacement_end], "A");
    assert!(mid.items.iter().any(|item| item.label == "A"));

    let after_keyword = complete_with_values("protein ", 8, &aliases, Some(&chains));
    assert!(after_keyword.items.is_empty());
}

#[test]
fn completion_value_candidates_borrow_without_repeating_registry_items() {
    let aliases = QueryAliases::new();
    let chains = Chains { labels: &["chain"] };
    let result = complete_with_values("chain ", 6, &aliases, Some(&chains));
    // `chain` is also a column label; the value candidate must not duplicate
    // it, and the structure-derived entry is typed as a value.
    let value_items = result
        .items
        .iter()
        .filter(|item| item.kind == CompletionKind::Value)
        .map(|item| item.label.as_str())
        .collect::<Vec<_>>();
    assert_eq!(value_items, ["chain"]);
}

#[test]
fn completion_at_fresh_token_offers_registry_not_values() {
    let aliases = QueryAliases::new();
    let chains = Chains { labels: &["A"] };
    let fresh = complete_with_values("", 0, &aliases, Some(&chains));
    assert!(
        fresh
            .items
            .iter()
            .all(|item| item.kind != CompletionKind::Value)
    );
    assert!(fresh.items.iter().any(|item| item.label == "protein"));
}
