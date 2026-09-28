//! Completion for the selection language.
//!
//! This module is the single implementation of completion: candidates come
//! from the language registry and the alias registry, and structure-derived
//! values are lent through [`StructureValues`] rather than copied into a
//! second vocabulary. Nothing here parses the language a second time.

use crate::api::QueryAliases;
use crate::lexer::lex;
use molframe_core::diagnostic::Diagnostic;
use std::collections::HashSet;

/// The category of a completion item.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CompletionKind {
    /// A language keyword such as `protein` or `and`.
    Keyword,
    /// A structure-property column such as `chain` or `bfactor`.
    Column,
    /// A built-in structural expansion such as `within` or `byres`.
    Macro,
    /// A named query supplied through the alias registry.
    Alias,
    /// A structure-derived value for the preceding column, such as a chain
    /// identifier after `chain`.
    Value,
}

/// One supported completion candidate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompletionItem {
    /// Text inserted when this item is accepted.
    pub label: String,
    /// Candidate category used by interactive clients for presentation.
    pub kind: CompletionKind,
}

/// Lent, structure-derived completion values.
///
/// Implementors expose identifiers from native indexes as borrowed strings;
/// completion never builds a second vocabulary or copies the structure's
/// columns. A value the structure does not carry is simply omitted.
pub trait StructureValues {
    /// Distinct chain identifiers, in native order. Both namespaces are
    /// offered when a structure stores them separately.
    fn chain_labels(&self) -> Box<dyn Iterator<Item = &str> + '_>;

    /// Distinct residue names, in native order.
    fn residue_names(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        Box::new(std::iter::empty())
    }

    /// Distinct atom names, in native order.
    fn atom_names(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        Box::new(std::iter::empty())
    }
}

impl StructureValues for molframe_core::structure::Structure {
    fn chain_labels(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        let mut seen = HashSet::new();
        Box::new(
            self.data()
                .chains()
                .flat_map(|chain| [chain.label(), chain.auth_label()])
                .flatten()
                .filter(move |label| seen.insert(*label)),
        )
    }

    fn residue_names(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        let mut seen = HashSet::new();
        Box::new(
            self.data()
                .residues()
                .filter_map(molframe_core::structure::ResidueRef::name)
                .filter(move |name| seen.insert(*name)),
        )
    }

    fn atom_names(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        let mut seen = HashSet::new();
        Box::new(
            self.data()
                .atoms()
                .filter_map(molframe_core::structure::AtomRef::name)
                .filter(move |name| seen.insert(*name)),
        )
    }
}

/// Completion response with a byte-addressed replacement range.
#[derive(Clone, Debug)]
pub struct CompletionResult {
    /// Start of the token to replace, in UTF-8 bytes.
    pub replacement_start: usize,
    /// End of the token to replace, in UTF-8 bytes.
    pub replacement_end: usize,
    /// Candidates in deterministic registry order.
    pub items: Vec<CompletionItem>,
    /// Diagnostics produced while lexing the prefix.
    pub diagnostics: Vec<Diagnostic>,
}

const KEYWORDS: &[&str] = &[
    "and",
    "or",
    "not",
    "within",
    "beyond",
    "around",
    "sphzone",
    "sphlayer",
    "isolayer",
    "cyzone",
    "cylayer",
    "point",
    "of",
    "byres",
    "same",
    "as",
    "group",
    "connected",
    "bonded",
    "by",
    "to",
];

const COLUMNS: &[&str] = &[
    "index",
    "bynum",
    "id",
    "resindex",
    "chainindex",
    "modelindex",
    "chain",
    "label_chain",
    "auth_chain",
    "resid",
    "label_resid",
    "auth_resid",
    "resname",
    "label_resname",
    "auth_resname",
    "name",
    "label_name",
    "auth_name",
    "altloc",
    "model",
    "entity",
    "entity_type",
    "element",
    "segid",
    "icode",
    "record_type",
    "x",
    "y",
    "z",
    "mass",
    "charge",
    "formalcharge",
    "radius",
    "bfactor",
    "occupancy",
    "plddt",
    "pae",
];

const MACROS: &[&str] = &[
    "protein",
    "backbone",
    "sidechain",
    "nucleic",
    "nucleicbackbone",
    "nucleicbase",
    "nucleicsugar",
    "water",
    "ion",
    "lipid",
    "saccharide",
    "hetero",
    "hydrogen",
    "heavy",
    "polymer",
    "ligand",
    "aromatic",
];

/// Columns whose operands are identifiers a structure can lend, mapped to the
/// [`StructureValues`] method that supplies them.
fn value_column_candidates<'a>(
    column: &str,
    values: &'a dyn StructureValues,
) -> Box<dyn Iterator<Item = &'a str> + 'a> {
    match column {
        "chain" | "label_chain" | "auth_chain" => values.chain_labels(),
        "resname" | "label_resname" | "auth_resname" => values.residue_names(),
        "name" | "label_name" | "auth_name" => values.atom_names(),
        _ => Box::new(std::iter::empty()),
    }
}

/// The identifier immediately before the token being completed, if any.
///
/// Delimiters that [`token_fragment`] treats as boundaries are skipped, so
/// `(chain ` still resolves to the `chain` column.
fn previous_word(before: &str) -> Option<&str> {
    let trimmed = before.trim_end_matches(|character: char| {
        character.is_whitespace() || matches!(character, '(' | ')' | ',')
    });
    let word_start = trimmed
        .char_indices()
        .rev()
        .find(|(_, character)| character.is_whitespace() || matches!(character, '(' | ')' | ','))
        .map_or(0, |(index, character)| index + character.len_utf8());
    let word = &trimmed[word_start..];
    (!word.is_empty()).then_some(word)
}

/// Computes completion using a UTF-8 byte cursor.
///
/// `cursor` is a byte offset, not a character index. A cursor in the middle
/// of a UTF-8 scalar is rejected rather than silently moving to a different
/// replacement range. The returned range always points into `source`.
///
/// After a value column such as `chain`, `resname`, or `name`, candidates
/// come exclusively from `values` when it is supplied; language keywords are
/// never mixed with structure identifiers. `$` references are completed
/// exclusively from `aliases`. At a fresh token after a finished term, no
/// candidates are offered until the caller types a character.
#[must_use]
pub fn complete(
    source: &str,
    cursor: usize,
    aliases: &QueryAliases,
    values: Option<&dyn StructureValues>,
) -> CompletionResult {
    let mut diagnostics = Vec::new();
    if cursor > source.len() || !source.is_char_boundary(cursor) {
        diagnostics.push(crate::lexer::syntax(
            crate::lexer::end_of(source),
            "completion cursor must be a UTF-8 byte boundary",
        ));
        // Keep the response range valid even for an interior UTF-8 offset.
        let boundary = source
            .char_indices()
            .map(|(index, _)| index)
            .take_while(|&index| index < cursor)
            .last()
            .map_or(0, |value| value);
        return CompletionResult {
            replacement_start: boundary,
            replacement_end: boundary,
            items: Vec::new(),
            diagnostics,
        };
    }

    let prefix = &source[..cursor];
    let (start, fragment) = token_fragment(prefix);
    let after_operator = prefix
        .trim_end()
        .chars()
        .next_back()
        .is_some_and(|character| matches!(character, '=' | '<' | '>'));
    // Dollar-prefixed names are a distinct namespace. Restrict the check to
    // the current token so an earlier dollar cannot change completion mode,
    // and compare aliases without the sigil.
    let dollar_start = prefix.rfind('$').filter(|&index| {
        prefix[index + 1..]
            .chars()
            .all(|character| character.is_ascii_alphanumeric() || character == '_')
            && (index == 0
                || prefix[..index]
                    .chars()
                    .next_back()
                    .is_some_and(|character| {
                        character.is_whitespace() || matches!(character, '(' | ')' | ',')
                    }))
    });
    let fresh_start = prefix.trim().is_empty();
    let mut items = Vec::new();
    if let Some(dollar) = dollar_start {
        let alias_fragment = &prefix[dollar + 1..];
        for name in aliases.names() {
            if name.starts_with(alias_fragment) {
                items.push(CompletionItem {
                    label: name.to_owned(),
                    kind: CompletionKind::Alias,
                });
            }
        }
    } else if !after_operator {
        let column = previous_word(&prefix[..start])
            .and_then(|word| crate::ast::Column::from_name(&word.to_ascii_lowercase()));
        if let (Some(column), Some(values)) = (column, values) {
            let mut seen = HashSet::new();
            for value in value_column_candidates(column_name_of(column), values) {
                if value.starts_with(fragment) && seen.insert(value) {
                    items.push(CompletionItem {
                        label: value.to_owned(),
                        kind: CompletionKind::Value,
                    });
                }
            }
        } else if !fragment.is_empty() || fresh_start {
            for &(candidate, kind) in &[
                (KEYWORDS, CompletionKind::Keyword),
                (COLUMNS, CompletionKind::Column),
                (MACROS, CompletionKind::Macro),
            ] {
                for &label in candidate {
                    if label.starts_with(fragment) {
                        items.push(CompletionItem {
                            label: label.to_owned(),
                            kind,
                        });
                    }
                }
            }
        }
    }
    if let Err(finding) = lex(prefix) {
        diagnostics.push(finding);
    }
    CompletionResult {
        replacement_start: start,
        replacement_end: cursor,
        items,
        diagnostics,
    }
}

/// The registry label of a parsed column.
///
/// `Column::from_name` already canonicalised the input, so matching the
/// canonical spelling is exact.
fn column_name_of(column: crate::ast::Column) -> &'static str {
    match column {
        crate::ast::Column::Chain => "chain",
        crate::ast::Column::LabelChain => "label_chain",
        crate::ast::Column::AuthChain => "auth_chain",
        crate::ast::Column::ResidueName => "resname",
        crate::ast::Column::LabelResidueName => "label_resname",
        crate::ast::Column::AuthResidueName => "auth_resname",
        crate::ast::Column::AtomName => "name",
        crate::ast::Column::LabelAtomName => "label_name",
        crate::ast::Column::AuthAtomName => "auth_name",
        _ => "",
    }
}

fn token_fragment(source: &str) -> (usize, &str) {
    let start = source
        .char_indices()
        .rev()
        .find_map(|(index, character)| {
            (character.is_whitespace() || matches!(character, '(' | ')' | ','))
                .then_some(index + character.len_utf8())
        })
        .map_or(0, |value| value);
    (start, &source[start..])
}

#[cfg(test)]
#[path = "completion_tests.rs"]
mod tests;
