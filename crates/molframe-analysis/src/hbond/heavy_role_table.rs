//! Donor and acceptor typing rules for heavy atoms, from a table file.
//!
//! A [`HeavyRoleTable`] says which residue names are protein, nucleic acid or
//! water, and what role each named atom of them plays. The typing code in
//! `heavy_roles` only asks it, and no rows live in code: the table is selected by
//! name like every other. The text format is the shared table format, headed by
//! the four identity lines of every table:
//!
//! ```text
//! name: lab-roles
//! version: 1
//! source: our curation of the ligand set
//! licence: CC0-1.0
//! amino: ALA GLY
//! nucleic: DA
//! water: HOH
//!
//! # residues   atoms   role      [condition]
//! @amino       N       donor
//! @amino       O       acceptor
//! DA           N6      donor
//! ```
//!
//! `amino:`, `nucleic:` and `water:` list the residue names of each class and
//! must come before a row that names the class as `@amino`, `@nucleic` or
//! `@water`. A row is a comma list of residues, a comma list of atoms (a trailing
//! `*` matches a prefix), a role (`none`, `donor`, `acceptor`, `both`) and an
//! optional condition: `chain-start` (only the first residue of a chain) or `weak`
//! (only when sulfur is admitted). Protonation variants are looked up under their
//! own name and then their parent's; for one residue and atom the first row whose
//! condition holds applies. An atom of a listed residue that no row names has no role;
//! a residue in no class falls through to the generic rule.
//!
//! Complexity: decoding is O(rows x atoms listed); a role lookup is O(1) by hash
//! plus O(prefix rows of the residue).

use crate::protonation::parent_residue;
use molframe_chem::{
    Provided, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};
use std::collections::{BTreeSet, HashMap};

pub(super) const DONOR: u8 = 1;
pub(super) const ACCEPTOR: u8 = 2;
pub(super) const BOTH: u8 = DONOR | ACCEPTOR;

const WHAT: &str = "heavy-atom role table";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Condition {
    Always,
    ChainStart,
    Weak,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Rule {
    condition: Condition,
    role: u8,
}

/// Donor and acceptor roles by residue and atom name.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HeavyRoleTable {
    amino: BTreeSet<String>,
    nucleic: BTreeSet<String>,
    water: BTreeSet<String>,
    exact: HashMap<String, HashMap<String, Vec<Rule>>>,
    prefix: HashMap<String, Vec<(String, Rule)>>,
}

impl TableKind for HeavyRoleTable {
    const KIND: &'static str = "heavy-roles";
    const EXTENSION: &'static str = "roles";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(WHAT, text_of(WHAT, bytes)?, &mut table)?;
        if table.exact.is_empty() && table.prefix.is_empty() {
            return Err(TableTextError::at(WHAT, 0, "the table has no row"));
        }
        Ok(Provided { identity, table })
    }
}

impl TableLines for HeavyRoleTable {
    fn header(&mut self, line: usize, key: &str, value: &str) -> Result<bool, TableTextError> {
        if !matches!(key, "amino" | "nucleic" | "water") {
            return Ok(false);
        }
        declare_class(self, line, key, value)?;
        Ok(true)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        read_row(self, line, text)
    }
}

impl HeavyRoleTable {
    /// The role bits of `atom` in `component`, or `None` when the residue is in no
    /// class, so the generic rule applies.
    pub(super) fn role(
        &self,
        component: &str,
        atom: &str,
        starts_chain: bool,
        weak: bool,
    ) -> Option<u8> {
        let parent = parent_residue(component);
        let class = if self.water.contains(component) || self.nucleic.contains(component) {
            true
        } else {
            self.amino.contains(parent)
        };
        if !class {
            return None;
        }
        let find = |residue: &str| self.rule(residue, atom, starts_chain, weak);
        Some(
            find(component)
                .or_else(|| (parent != component).then(|| find(parent)).flatten())
                .unwrap_or(0),
        )
    }

    fn rule(&self, residue: &str, atom: &str, starts_chain: bool, weak: bool) -> Option<u8> {
        let holds = |rule: &&Rule| match rule.condition {
            Condition::Always => true,
            Condition::ChainStart => starts_chain,
            Condition::Weak => weak,
        };
        let exact = self
            .exact
            .get(residue)
            .and_then(|atoms| atoms.get(atom))
            .and_then(|rules| rules.iter().find(holds));
        if let Some(rule) = exact {
            return Some(rule.role);
        }
        self.prefix
            .get(residue)?
            .iter()
            .filter(|(prefix, _)| atom.starts_with(prefix.as_str()))
            .map(|(_, rule)| rule)
            .find(holds)
            .map(|rule| rule.role)
    }
}

fn declare_class(
    table: &mut HeavyRoleTable,
    line: usize,
    key: &str,
    value: &str,
) -> Result<(), TableTextError> {
    let names = value
        .split([',', ' '])
        .filter(|name| !name.is_empty())
        .map(str::to_owned);
    let class = match key {
        "amino" => &mut table.amino,
        "nucleic" => &mut table.nucleic,
        "water" => &mut table.water,
        other => {
            return Err(TableTextError::at(
                WHAT,
                line,
                format!("unknown header `{other}`"),
            ));
        }
    };
    if !class.is_empty() {
        return Err(TableTextError::at(
            WHAT,
            line,
            format!("`{key}` is declared twice"),
        ));
    }
    class.extend(names);
    if class.is_empty() {
        return Err(TableTextError::at(
            WHAT,
            line,
            format!("`{key}` has no value"),
        ));
    }
    Ok(())
}

fn read_row(table: &mut HeavyRoleTable, line: usize, text: &str) -> Result<(), TableTextError> {
    let words: Vec<&str> = text.split_whitespace().collect();
    let [residues, atoms, role, condition @ ..] = words.as_slice() else {
        return Err(TableTextError::at(
            WHAT,
            line,
            "expected `residues atoms role [condition]`",
        ));
    };
    let role = match *role {
        "none" => 0,
        "donor" => DONOR,
        "acceptor" => ACCEPTOR,
        "both" => BOTH,
        other => {
            return Err(TableTextError::at(
                WHAT,
                line,
                format!("unknown role `{other}`; expected none, donor, acceptor or both"),
            ));
        }
    };
    let condition = match condition {
        [] => Condition::Always,
        ["chain-start"] => Condition::ChainStart,
        ["weak"] => Condition::Weak,
        _ => {
            return Err(TableTextError::at(
                WHAT,
                line,
                "the condition must be `chain-start` or `weak`",
            ));
        }
    };
    let template = Rule { condition, role };
    let mut residue_names = Vec::new();
    for token in residues.split(',').filter(|token| !token.is_empty()) {
        residue_names.extend(expand(table, line, token)?);
    }
    for residue in &residue_names {
        for atom in atoms.split(',').filter(|atom| !atom.is_empty()) {
            if let Some(stem) = atom.strip_suffix('*') {
                table
                    .prefix
                    .entry(residue.clone())
                    .or_default()
                    .push((stem.to_owned(), template));
            } else {
                table
                    .exact
                    .entry(residue.clone())
                    .or_default()
                    .entry(atom.to_owned())
                    .or_default()
                    .push(template);
            }
        }
    }
    Ok(())
}

fn expand(table: &HeavyRoleTable, line: usize, token: &str) -> Result<Vec<String>, TableTextError> {
    let class = match token {
        "@amino" => &table.amino,
        "@nucleic" => &table.nucleic,
        "@water" => &table.water,
        name if name.starts_with('@') => {
            return Err(TableTextError::at(
                WHAT,
                line,
                format!("unknown class `{name}`; expected @amino, @nucleic or @water"),
            ));
        }
        name => return Ok(vec![name.to_owned()]),
    };
    if class.is_empty() {
        return Err(TableTextError::at(
            WHAT,
            line,
            format!("`{token}` is used before its class is declared"),
        ));
    }
    Ok(class.iter().cloned().collect())
}

#[cfg(test)]
#[path = "heavy_role_table_tests.rs"]
mod tests;
