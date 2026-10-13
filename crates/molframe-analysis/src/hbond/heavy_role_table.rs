//! Donor and acceptor typing rules for heavy atoms, built-in or the caller's own.
//!
//! A [`HeavyRoleTable`] says which residue names are protein, nucleic acid or
//! water, and what role each named atom of them plays. The typing code in
//! `heavy_roles` only asks it; the default rows live in `heavy_role_defaults`.
//! A caller's table has the same text format, headed by the four identity lines
//! of every table:
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

use crate::hbond::heavy_role_defaults::{HEAVY_ROLES_DEFAULT, HEAVY_ROLES_DEFAULT_TEXT};
use crate::protonation::parent_residue;
use crate::table_text::{IdentityFields, TableTextError};
use molframe_chem::TableIdentity;
use std::collections::{BTreeSet, HashMap};
use std::sync::OnceLock;

pub(super) const DONOR: u8 = 1;
pub(super) const ACCEPTOR: u8 = 2;
pub(super) const BOTH: u8 = DONOR | ACCEPTOR;

const WHAT: &str = "heavy-atom role table";

/// Why a role table was refused.
pub type HeavyRoleTableError = TableTextError;

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

/// Donor and acceptor roles by residue and atom name, with the identity of their source.
#[derive(Clone, Debug, PartialEq)]
pub struct HeavyRoleTable {
    identity: TableIdentity,
    amino: BTreeSet<String>,
    nucleic: BTreeSet<String>,
    water: BTreeSet<String>,
    exact: HashMap<String, HashMap<String, Vec<Rule>>>,
    prefix: HashMap<String, Vec<(String, Rule)>>,
}

impl HeavyRoleTable {
    /// The built-in table, parsed once and borrowed for the life of the process.
    #[must_use]
    pub fn builtin() -> &'static Self {
        static BUILTIN: OnceLock<HeavyRoleTable> = OnceLock::new();
        BUILTIN.get_or_init(|| {
            // The default text is pinned by a test, so a failure here cannot reach a
            // caller; an empty table is the safe reading of one.
            build(HEAVY_ROLES_DEFAULT_TEXT, Some(HEAVY_ROLES_DEFAULT))
                .unwrap_or_else(|_| Self::empty(HEAVY_ROLES_DEFAULT))
        })
    }

    fn empty(identity: TableIdentity) -> Self {
        Self {
            identity,
            amino: BTreeSet::new(),
            nucleic: BTreeSet::new(),
            water: BTreeSet::new(),
            exact: HashMap::new(),
            prefix: HashMap::new(),
        }
    }

    /// Parses and validates a table.
    ///
    /// # Errors
    ///
    /// [`HeavyRoleTableError`] for a missing, unknown or repeated header, an unknown
    /// role or condition, a class used before it is declared, or a malformed row.
    pub fn parse(text: &str) -> Result<Self, HeavyRoleTableError> {
        build(text, None)
    }

    /// The declared identity of this table.
    #[must_use]
    pub const fn identity(&self) -> &TableIdentity {
        &self.identity
    }

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

fn build(text: &str, fixed: Option<TableIdentity>) -> Result<HeavyRoleTable, TableTextError> {
    let mut fields = IdentityFields::default();
    let mut table = HeavyRoleTable::empty(HEAVY_ROLES_DEFAULT);
    let mut rows = 0_usize;
    for (index, raw) in text.lines().enumerate() {
        let number = index + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some((key, value)) = line.split_once(':') {
            let (key, value) = (key.trim(), value.trim());
            if !fields.set(WHAT, number, key, value)? {
                declare_class(&mut table, number, key, value)?;
            }
        } else {
            read_row(&mut table, number, line)?;
            rows += 1;
        }
    }
    if rows == 0 {
        return Err(TableTextError::at(WHAT, 0, "the table has no row"));
    }
    table.identity = match fixed {
        Some(identity) => identity,
        None => fields.finish(WHAT)?,
    };
    Ok(table)
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
                format!(
                    "unknown header `{other}`; expected name, version, source, licence, amino, nucleic or water"
                ),
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
