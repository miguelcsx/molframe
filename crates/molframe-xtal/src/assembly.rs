//! Normalised biological-assembly definitions.

use crate::OperExpression;
use molframe_core::{Code, Diagnostic};
use molframe_geom::Rigid;
use std::collections::BTreeMap;

/// Stable extension key used to attach assembly metadata to a structure.
pub const ASSEMBLIES_EXTENSION: &str = "molframe.xtal.assemblies.v1";

/// One Cartesian transform named by the source dictionary.
#[derive(Clone, Debug, PartialEq)]
pub struct Operator {
    /// Source operator identifier.
    pub id: Box<str>,
    /// Rotation followed by translation in Cartesian coordinates.
    pub transform: Rigid,
}

impl Operator {
    /// An operator from a row-major rotation and a translation, both Cartesian.
    #[must_use]
    pub fn from_matrix(
        id: impl Into<Box<str>>,
        rotation: [[f64; 3]; 3],
        translation: [f64; 3],
    ) -> Self {
        Self {
            id: id.into(),
            transform: Rigid::new(rotation, translation),
        }
    }
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

    /// Builds a set from definitions and the operators they name, for sources
    /// other than mmCIF (such as a legacy PDB `REMARK 350`).
    ///
    /// # Errors
    ///
    /// Returns [`Code::E6012`] for a repeated assembly or operator identifier,
    /// and [`Code::E6013`] for a generator naming an operator that was not
    /// supplied.
    pub fn from_definitions(
        assemblies: impl IntoIterator<Item = AssemblyDef>,
        operators: impl IntoIterator<Item = Operator>,
    ) -> Result<Self, Diagnostic> {
        let mut set = Self::default();
        for operator in operators {
            let id = operator.id.clone();
            if !crate::lower::valid_rotation(&operator.transform)
                || set.operators.insert(id.clone(), operator).is_some()
            {
                return Err(Diagnostic::new(Code::E6012).with_context("operator", id));
            }
        }
        for assembly in assemblies {
            let id = assembly.id.clone();
            let named = assembly
                .generators
                .iter()
                .flat_map(|generator| generator.oper_expression.factors())
                .flat_map(|factor| factor.iter());
            for operator in named {
                if !set.operators.contains_key(operator) {
                    return Err(Diagnostic::new(Code::E6013)
                        .with_context("assembly", id)
                        .with_context("operator", operator.to_string()));
                }
            }
            if set.assemblies.insert(id.clone(), assembly).is_some() {
                return Err(Diagnostic::new(Code::E6012).with_context("assembly", id));
            }
        }
        Ok(set)
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
