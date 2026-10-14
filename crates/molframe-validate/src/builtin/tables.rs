//! Bond, angle and chi-atom tables.
//!
//! Both are tab-separated text under the shared table format and are decoded
//! from any source by the same code. A malformed row is an error, never skipped.
//! Rows are `bond|angle|link residue atoms mean sd [source]` for geometry and
//! `residue table chi | chi` for chi atoms.
//!
//! Complexity: decoding is O(rows) once, then shared through an `Arc`; a
//! residue lookup is O(1) by hash and a link lookup is O(link terms).

use std::collections::HashMap;

use molframe_chem::{
    Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};

/// A bond length or bond angle target with its standard deviation.
#[derive(Clone, Debug, PartialEq)]
pub struct Target {
    /// Atom names in bonded order; two for a bond, three for an angle.
    pub names: Vec<String>,
    /// Ideal value, in ångström or degrees.
    pub mean: f64,
    /// Standard deviation of the ideal value.
    pub sigma: f64,
}

impl Target {
    /// Whether this target is an angle rather than a bond length.
    #[must_use]
    pub fn is_angle(&self) -> bool {
        self.names.len() == 3
    }
}

/// Targets keyed by residue, plus the peptide-link terms.
#[derive(Debug, Default)]
pub struct GeometryTargets {
    residue: HashMap<String, Vec<Target>>,
    link: Vec<Target>,
    link_special: HashMap<String, Vec<Target>>,
}

impl GeometryTargets {
    /// Targets inside one residue, by its parent residue name.
    #[must_use]
    pub fn residue(&self, name: &str) -> &[Target] {
        match self.residue.get(name) {
            Some(targets) => targets,
            None => &[],
        }
    }

    /// Peptide-link terms for a residue that follows a carbonyl, with the
    /// entries keyed by that residue replacing the generic ones of the same atoms.
    pub fn link<'a>(&'a self, following: &str) -> impl Iterator<Item = &'a Target> + 'a {
        let special = self.link_special.get(following);
        self.link.iter().map(move |generic| {
            let replacement = special.and_then(|terms| {
                terms
                    .iter()
                    .find(|candidate| candidate.names == generic.names)
            });
            match replacement {
                Some(term) => term,
                None => generic,
            }
        })
    }

    /// Total number of targets.
    #[must_use]
    pub fn len(&self) -> usize {
        self.residue.values().map(Vec::len).sum::<usize>()
            + self.link.len()
            + self.link_special.values().map(Vec::len).sum::<usize>()
    }

    /// Whether no targets were loaded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

const GEOMETRY: &str = "geometry table";
const CHI: &str = "chi table";

impl TableLines for GeometryTargets {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let fields: Vec<&str> = text.split('\t').collect();
        let (&[kind, residue, atoms, mean, sigma] | &[kind, residue, atoms, mean, sigma, _]) =
            fields.as_slice()
        else {
            return Err(TableTextError::at(
                GEOMETRY,
                line,
                "expected `kind residue atoms mean sd [source]` separated by tabs",
            ));
        };
        let names: Vec<String> = atoms.split('-').map(str::to_owned).collect();
        if names.iter().any(String::is_empty) || !matches!(names.len(), 2 | 3) {
            return Err(TableTextError::at(
                GEOMETRY,
                line,
                format!("`{atoms}` is not two or three atom names joined by `-`"),
            ));
        }
        let mean = Row::new(GEOMETRY, line).finite::<f64>(mean)?;
        let sigma = Row::new(GEOMETRY, line).finite::<f64>(sigma)?;
        if sigma <= 0.0 {
            return Err(TableTextError::at(GEOMETRY, line, "sd must be positive"));
        }
        let target = Target { names, mean, sigma };
        match (kind, residue) {
            ("bond", name) if target.names.len() == 2 && name != "*" => {
                self.residue
                    .entry(name.to_owned())
                    .or_default()
                    .push(target);
            }
            ("angle", name) if target.names.len() == 3 && name != "*" => {
                self.residue
                    .entry(name.to_owned())
                    .or_default()
                    .push(target);
            }
            ("link", "*") => self.link.push(target),
            ("link", name) => self
                .link_special
                .entry(name.to_owned())
                .or_default()
                .push(target),
            _ => {
                return Err(TableTextError::at(
                    GEOMETRY,
                    line,
                    format!("`{kind}` with residue `{residue}` and `{atoms}` is not a valid row"),
                ));
            }
        }
        Ok(())
    }
}

impl TableKind for GeometryTargets {
    const KIND: &'static str = "geometry";
    const EXTENSION: &'static str = "geometry";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(GEOMETRY, text_of(GEOMETRY, bytes)?, &mut table)?;
        if table.is_empty() {
            return Err(TableTextError::at(GEOMETRY, 0, "the table has no rows"));
        }
        Ok(Provided { identity, table })
    }
}

/// The atoms that define each chi angle of one residue type.
#[derive(Clone, Debug, PartialEq)]
pub struct ChiDefinition {
    /// Name of the rotamer table to score against.
    pub table: String,
    /// Four atom names per chi, in order chi1, chi2, ...
    pub chi: Vec<[String; 4]>,
}

/// Chi definitions keyed by parent residue name.
#[derive(Clone, Debug, Default)]
pub struct ChiDefinitions(HashMap<String, ChiDefinition>);

impl ChiDefinitions {
    /// The chi definition of one residue type, by parent name.
    #[must_use]
    pub fn get(&self, residue: &str) -> Option<&ChiDefinition> {
        self.0.get(residue)
    }

    /// Number of residue types defined.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Whether nothing is defined.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl TableLines for ChiDefinitions {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let fields: Vec<&str> = text.split('\t').collect();
        let &[residue, table, angles] = fields.as_slice() else {
            return Err(TableTextError::at(
                CHI,
                line,
                "expected `residue table chi | chi` separated by tabs",
            ));
        };
        let mut chi = Vec::new();
        for part in angles.split('|') {
            let names: Vec<String> = part.split_whitespace().map(str::to_owned).collect();
            let Ok(four) = <[String; 4]>::try_from(names) else {
                return Err(TableTextError::at(
                    CHI,
                    line,
                    format!("`{}` is not four atom names", part.trim()),
                ));
            };
            chi.push(four);
        }
        let definition = ChiDefinition {
            table: table.to_owned(),
            chi,
        };
        if self.0.insert(residue.to_owned(), definition).is_some() {
            return Err(TableTextError::at(
                CHI,
                line,
                format!("residue `{residue}` is defined twice"),
            ));
        }
        Ok(())
    }
}

impl TableKind for ChiDefinitions {
    const KIND: &'static str = "chi";
    const EXTENSION: &'static str = "chi";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(CHI, text_of(CHI, bytes)?, &mut table)?;
        if table.is_empty() {
            return Err(TableTextError::at(CHI, 0, "the table has no rows"));
        }
        Ok(Provided { identity, table })
    }
}

#[cfg(test)]
#[path = "tables_tests.rs"]
mod tests;
