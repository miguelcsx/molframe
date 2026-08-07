//! Normalised biological-assembly definitions.

use crate::OperExpression;
use pdbiox_geom::Rigid;
use std::collections::BTreeMap;

/// Stable extension key used to attach assembly metadata to a structure.
pub const ASSEMBLIES_EXTENSION: &str = "pdbiox.xtal.assemblies.v1";

/// One Cartesian transform named by the source dictionary.
#[derive(Clone, Debug, PartialEq)]
pub struct Operator {
    /// Source operator identifier.
    pub id: Box<str>,
    /// Rotation followed by translation in Cartesian coordinates.
    pub transform: Rigid,
}

/// One rule applying an operator expression to label-asym identifiers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Generator {
    /// Ordered operator product.
    pub oper_expression: OperExpression,
    /// `_atom_site.label_asym_id` values, in declared order.
    pub asym_ids: Box<[Box<str>]>,
}

/// One biological assembly declared by the entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AssemblyDef {
    /// Source assembly identifier.
    pub id: Box<str>,
    /// Depositor or software description.
    pub details: Option<Box<str>>,
    /// Method used to determine the assembly.
    pub method: Option<Box<str>>,
    /// Declared oligomeric count.
    pub oligomeric: Option<u32>,
    /// Rules that generate its chain instances.
    pub generators: Vec<Generator>,
}

/// Assembly definitions and their shared operator dictionary.
#[derive(Clone, Debug, Default)]
pub struct AssemblySet {
    assemblies: BTreeMap<Box<str>, AssemblyDef>,
    operators: BTreeMap<Box<str>, Operator>,
}

impl AssemblySet {
    pub(crate) fn from_parts(
        assemblies: BTreeMap<Box<str>, AssemblyDef>,
        operators: BTreeMap<Box<str>, Operator>,
    ) -> Self {
        Self {
            assemblies,
            operators,
        }
    }

    /// One assembly by its source identifier.
    #[must_use]
    pub fn get(&self, id: &str) -> Option<&AssemblyDef> {
        self.assemblies.get(id)
    }

    /// One operator by its source identifier.
    #[must_use]
    pub fn operator(&self, id: &str) -> Option<&Operator> {
        self.operators.get(id)
    }

    /// Assemblies in lexical identifier order.
    pub fn assemblies(&self) -> impl Iterator<Item = &AssemblyDef> {
        self.assemblies.values()
    }

    /// Operators in lexical identifier order.
    pub fn operators(&self) -> impl Iterator<Item = &Operator> {
        self.operators.values()
    }

    /// Number of declared assemblies.
    #[must_use]
    pub fn len(&self) -> usize {
        self.assemblies.len()
    }

    /// Whether no assembly is declared.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.assemblies.is_empty()
    }
}
