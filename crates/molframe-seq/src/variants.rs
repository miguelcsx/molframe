//! Tables of amino-acid substitutions with a score each (a deep mutational scan's measured
//! effects, a language model's predictions, a stability calculator's ΔΔG), and the position
//! summaries and numbering checks built on them in [`crate::variant_effects`].

use crate::residue_letters::one_letter;
use thiserror::Error;

/// Header words that name the column of substitutions.
const MUTANT_COLUMNS: [&str; 7] = [
    "mutant",
    "variant",
    "mutation",
    "mutations",
    "substitution",
    "substitutions",
    "aa_change",
];

/// One substitution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Substitution {
    /// The residue in the wild type.
    pub wild_type: char,
    /// Its position in the table's numbering.
    pub position: i32,
    /// The residue it became (`*` for a stop).
    pub mutant: char,
}

/// One variant: one or several substitutions, and the values the table gives it.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Variant {
    /// As written.
    pub label: String,
    /// The substitutions.
    pub changes: Vec<Substitution>,
    /// One value per column of the table (`None` where the cell was empty).
    pub values: Vec<Option<f64>>,
}

/// What a variant table said.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct VariantTable {
    /// The file's name.
    pub source: String,
    /// The numeric columns, in file order.
    pub columns: Vec<String>,
    /// The variants.
    pub variants: Vec<Variant>,
    /// Rows left out, and why.
    pub notes: Vec<String>,
}

/// Why a variant table could not be read. The messages are complete sentences.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum VariantTableError {
    /// The text has no header line.
    #[error("the file is empty")]
    Empty,
    /// No header names the substitutions.
    #[error("no column of substitutions (mutant, variant, mutation)")]
    NoSubstitutionColumn,
    /// No column holds only numbers.
    #[error("no numeric column: a variant needs a score")]
    NoNumericColumn,
    /// No row holds a readable substitution.
    #[error("no row holds a substitution like A123G")]
    NoReadableRow,
}

fn three_to_one(name: &str) -> Option<char> {
    let letter = one_letter(name);
    (letter != 'X').then_some(letter)
}

/// Reads one substitution: `A123G`, `Ala123Gly`, `p.Ala123Gly`, `A123*`.
#[must_use]
pub fn parse_substitution(text: &str) -> Option<Substitution> {
    let trimmed = text.trim().trim_start_matches("p.");
    let digits_at = trimmed.find(|c: char| c.is_ascii_digit())?;
    let (head, rest) = trimmed.split_at(digits_at);
    let end = match rest.find(|c: char| !c.is_ascii_digit()) {
        Some(end) => end,
        None => rest.len(),
    };
    let (number, tail) = rest.split_at(end);
    let position: i32 = number.parse().ok()?;
    let letter = |name: &str| -> Option<char> {
        match name.len() {
            1 => {
                let c = name.chars().next()?.to_ascii_uppercase();
                (c.is_ascii_alphabetic() || c == '*').then_some(c)
            }
            3 => three_to_one(name),
            _ => None,
        }
    };
    let tail = if tail.eq_ignore_ascii_case("ter") {
        "*"
    } else {
        tail
    };
    Some(Substitution {
        wild_type: letter(head)?,
        position,
        mutant: letter(tail)?,
    })
}

/// Splits one delimited row into trimmed cells. A double-quoted cell may hold the delimiter, and
/// a doubled quote inside one is a literal quote.
#[must_use]
pub fn split_delimited_row(line: &str, delimiter: char) -> Vec<String> {
    let mut cells = Vec::new();
    let mut cell = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' if quoted && chars.peek() == Some(&'"') => {
                cell.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            c if c == delimiter && !quoted => cells.push(std::mem::take(&mut cell)),
            c => cell.push(c),
        }
    }
    cells.push(cell);
    cells.into_iter().map(|c| c.trim().to_owned()).collect()
}

/// Reads a table of variants, comma- or tab-delimited.
///
/// # Errors
/// [`VariantTableError`] when the text has no header, no column of substitutions, no numeric
/// column, or no readable variant.
pub fn parse_variants(text: &str, source: &str) -> Result<VariantTable, VariantTableError> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines.next().ok_or(VariantTableError::Empty)?;
    let delimiter = if header.matches('\t').count() > header.matches(',').count() {
        '\t'
    } else {
        ','
    };
    let columns = split_delimited_row(header, delimiter);
    let mutant_at = columns
        .iter()
        .position(|c| MUTANT_COLUMNS.contains(&c.to_ascii_lowercase().as_str()))
        .ok_or(VariantTableError::NoSubstitutionColumn)?;
    let body: Vec<Vec<String>> = lines.map(|l| split_delimited_row(l, delimiter)).collect();
    let mut numeric = Vec::new();
    for (k, name) in columns.iter().enumerate() {
        if k == mutant_at || name.is_empty() {
            continue;
        }
        let cells: Vec<&str> = body
            .iter()
            .filter_map(|r| r.get(k).map(String::as_str))
            .filter(|c| !c.is_empty())
            .collect();
        if !cells.is_empty()
            && cells
                .iter()
                .all(|c| c.parse::<f64>().is_ok_and(f64::is_finite))
        {
            numeric.push(k);
        }
    }
    if numeric.is_empty() {
        return Err(VariantTableError::NoNumericColumn);
    }
    let mut variants = Vec::new();
    let mut notes = Vec::new();
    for (n, row) in body.iter().enumerate() {
        let Some(label) = row.get(mutant_at).filter(|c| !c.is_empty()) else {
            continue;
        };
        let changes: Option<Vec<Substitution>> = label
            .split([':', ',', ';', '+'])
            .map(parse_substitution)
            .collect();
        let Some(changes) = changes else {
            if notes.len() < 20 {
                notes.push(format!(
                    "row {}: {label:?} is not a substitution like A123G",
                    n + 2
                ));
            }
            continue;
        };
        variants.push(Variant {
            label: label.clone(),
            changes,
            values: numeric
                .iter()
                .map(|k| row.get(*k).and_then(|c| c.parse::<f64>().ok()))
                .collect(),
        });
    }
    if variants.is_empty() {
        return Err(VariantTableError::NoReadableRow);
    }
    let skipped = body.len().saturating_sub(variants.len());
    if skipped > notes.len() {
        notes.push(format!("{skipped} rows in all were not read"));
    }
    Ok(VariantTable {
        source: source.to_owned(),
        columns: numeric
            .iter()
            .filter_map(|k| columns.get(*k).cloned())
            .collect(),
        variants,
        notes,
    })
}

#[cfg(test)]
#[path = "variants_tests.rs"]
mod tests;
