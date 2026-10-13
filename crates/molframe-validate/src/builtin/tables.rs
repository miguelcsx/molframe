//! Parsers for bond, angle and chi-atom tables.
//!
//! Both formats are tab-separated text, documented in `data/README.md`, and
//! are parsed from any source: the bundled files and user-supplied ones go
//! through the same code. Malformed lines are skipped rather than trusted.

use std::collections::HashMap;

/// Atom-name patterns that join one residue to the previous one.
const LINK_PATTERNS: [&str; 4] = ["C-N", "CA-C-N", "C-N-CA", "O-C-N"];

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
    proline_link: Vec<Target>,
}

impl GeometryTargets {
    /// Parses a bond and angle target table.
    #[must_use]
    pub fn parse(text: &str) -> Self {
        parse_geometry(text)
    }

    /// Targets inside one residue, by its parent residue name.
    #[must_use]
    pub fn residue(&self, name: &str) -> &[Target] {
        match self.residue.get(name) {
            Some(targets) => targets,
            None => &[],
        }
    }

    /// Peptide-link terms for a residue that follows a carbonyl, with the
    /// proline-specific entries replacing the generic ones of the same atoms.
    pub fn link<'a>(&'a self, following: &str) -> impl Iterator<Item = &'a Target> + 'a {
        let proline = following == "PRO";
        self.link.iter().map(move |generic| {
            if proline {
                match self
                    .proline_link
                    .iter()
                    .find(|special| special.names == generic.names)
                {
                    Some(special) => special,
                    None => generic,
                }
            } else {
                generic
            }
        })
    }

    /// Total number of targets.
    #[must_use]
    pub fn len(&self) -> usize {
        self.residue.values().map(Vec::len).sum::<usize>()
            + self.link.len()
            + self.proline_link.len()
    }

    /// Whether no targets were loaded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

fn parse_geometry(text: &str) -> GeometryTargets {
    let mut targets = GeometryTargets::default();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let [_, residue, atoms, mean, sigma, ..] = fields.as_slice() else {
            continue;
        };
        let (Ok(mean), Ok(sigma)) = (mean.parse::<f64>(), sigma.parse::<f64>()) else {
            continue;
        };
        let target = Target {
            names: atoms.split('-').map(str::to_owned).collect(),
            mean,
            sigma,
        };
        let link = LINK_PATTERNS.contains(atoms);
        match (*residue, link) {
            ("*", _) => targets.link.push(target),
            ("PRO", true) => targets.proline_link.push(target),
            (name, false) => targets
                .residue
                .entry(name.to_owned())
                .or_default()
                .push(target),
            _ => {}
        }
    }
    targets
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
    /// Parses a chi-atom table.
    #[must_use]
    pub fn parse(text: &str) -> Self {
        Self(parse_chi(text))
    }

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

fn parse_chi(text: &str) -> HashMap<String, ChiDefinition> {
    let mut map = HashMap::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut head = line.split('\t');
        let (Some(residue), Some(table), Some(rest)) = (head.next(), head.next(), head.next())
        else {
            continue;
        };
        let chi = rest
            .split('|')
            .filter_map(|part| {
                let names: Vec<String> = part.split_whitespace().map(str::to_owned).collect();
                <[String; 4]>::try_from(names).ok()
            })
            .collect();
        map.insert(
            residue.to_owned(),
            ChiDefinition {
                table: table.to_owned(),
                chi,
            },
        );
    }
    map
}

#[cfg(test)]
#[path = "tables_tests.rs"]
mod tests;
