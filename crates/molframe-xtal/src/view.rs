//! Lazy biological-assembly instances over an immutable structure snapshot.

use crate::numeric::usize_to_u32;
use crate::{ASSEMBLIES_EXTENSION, AssemblySet, DEFAULT_INSTANCE_LIMIT};
use molframe_core::{AtomIndex, ChainIndex, Code, Diagnostic, InstanceId, ModelIndex, Structure};
use molframe_geom::Rigid;
use std::ops::Range;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub(crate) struct InstanceRecord {
    pub(crate) source_chain: ChainIndex,
    pub(crate) transform: usize,
    pub(crate) instance_id: InstanceId,
    pub(crate) atoms: Range<u32>,
}

/// One generated chain without copied topology or coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChainInstance {
    /// Chain in the deposited asymmetric unit.
    pub source_chain: ChainIndex,
    /// Cartesian transform applied to this copy.
    pub transform: Rigid,
    /// Stable identifier of this generated chain copy.
    pub instance_id: InstanceId,
}

/// One generated atom without copied atom data or coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomInstance {
    /// Atom in the deposited asymmetric unit.
    pub source_atom: AtomIndex,
    /// Cartesian transform applied to this copy.
    pub transform: Rigid,
    /// Stable identifier shared with its generated chain copy.
    pub instance_id: InstanceId,
}

impl AtomInstance {
    /// Applies this instance's transform to a source position.
    #[must_use]
    pub fn apply(self, source_position: [f32; 3]) -> [f32; 3] {
        self.transform.apply(source_position)
    }
}

/// A biological assembly represented as chain-transform pairs.
#[derive(Clone, Debug)]
pub struct AssemblyView {
    pub(crate) structure: Structure,
    assembly_id: Box<str>,
    pub(crate) instances: Arc<[InstanceRecord]>,
    pub(crate) transforms: Arc<[Rigid]>,
}

impl AssemblyView {
    /// Builds one named assembly under the default instance limit.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an absent assembly, unknown operator or chain,
    /// or a generated chain count above the configured ceiling.
    pub fn new(
        structure: &Structure,
        assemblies: &AssemblySet,
        id: &str,
    ) -> Result<Self, Diagnostic> {
        Self::with_limit(structure, assemblies, id, DEFAULT_INSTANCE_LIMIT)
    }

    /// Builds one named assembly under an explicit chain-instance limit.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for unresolved references or excessive expansion.
    pub fn with_limit(
        structure: &Structure,
        assemblies: &AssemblySet,
        id: &str,
        limit: usize,
    ) -> Result<Self, Diagnostic> {
        let definition = assemblies
            .get(id)
            .ok_or_else(|| Diagnostic::new(Code::E6002).with_context("assembly", id))?;
        let mut instances = Vec::new();
        let mut transforms = Vec::new();
        for generator in &definition.generators {
            let chains = resolve_chains(structure, &generator.asym_ids)?;
            let generated = generator
                .oper_expression
                .combination_count()
                .checked_mul(chains.len())
                .and_then(|count| instances.len().checked_add(count))
                .filter(|count| *count <= limit)
                .ok_or_else(|| instance_limit(id, limit))?;
            let additional = generated
                .checked_sub(instances.len())
                .ok_or_else(|| instance_limit(id, limit))?;
            instances
                .try_reserve(additional)
                .map_err(|_| instance_limit(id, limit))?;
            transforms
                .try_reserve(additional)
                .map_err(|_| instance_limit(id, limit))?;
            generator
                .oper_expression
                .for_each_combination(|identifiers| {
                    if instances.len() >= generated {
                        return;
                    }
                    let transform = compose(assemblies, identifiers);
                    let Ok(transform) = transform else {
                        return;
                    };
                    let transform_index = transforms.len();
                    transforms.push(transform);
                    for (chain, atoms) in &chains {
                        instances.push(InstanceRecord {
                            source_chain: *chain,
                            transform: transform_index,
                            instance_id: InstanceId::new(usize_to_u32(instances.len())),
                            atoms: atoms.clone(),
                        });
                    }
                });
            if instances.len() != generated {
                return Err(first_unknown_operator(
                    assemblies,
                    &generator.oper_expression,
                ));
            }
        }
        Ok(Self {
            structure: structure.clone(),
            assembly_id: definition.id.clone(),
            instances: instances.into(),
            transforms: transforms.into(),
        })
    }

    /// Source assembly identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.assembly_id
    }

    /// Deposited structure snapshot from which instances are generated.
    #[must_use]
    pub const fn source(&self) -> &Structure {
        &self.structure
    }

    /// Number of generated chain copies.
    #[must_use]
    pub fn instance_count(&self) -> usize {
        self.instances.len()
    }

    /// Generated chains in assembly declaration order.
    pub fn chains(&self) -> impl Iterator<Item = ChainInstance> + '_ {
        self.instances.iter().map(|instance| ChainInstance {
            source_chain: instance.source_chain,
            transform: self.transforms[instance.transform],
            instance_id: instance.instance_id,
        })
    }

    /// Generated atoms in chain-instance then source-atom order.
    pub fn atoms(&self) -> impl Iterator<Item = AtomInstance> + '_ {
        self.instances.iter().flat_map(|instance| {
            let transform = self.transforms[instance.transform];
            instance.atoms.clone().map(move |atom| AtomInstance {
                source_atom: AtomIndex::new(atom),
                transform,
                instance_id: instance.instance_id,
            })
        })
    }

    /// Generated atom positions for one source model, preserving missingness.
    pub fn positions(
        &self,
        model: ModelIndex,
    ) -> impl Iterator<Item = (AtomInstance, Option<[f32; 3]>)> + '_ {
        let positions = self.structure.model_positions(model);
        self.atoms().map(move |instance| {
            let position = positions
                .and_then(|values| values.get(instance.source_atom.as_usize()))
                .copied()
                .map(|position| instance.apply(position));
            (instance, position)
        })
    }
}

/// Assembly access provided by the typed structure extension.
pub trait AssemblyExt {
    /// Attached assembly definitions, when the reader found any.
    fn assembly_set(&self) -> Option<&AssemblySet>;

    /// Creates a lazy view of one named biological assembly.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic if metadata is absent or references do not resolve.
    fn assembly(&self, id: &str) -> Result<AssemblyView, Diagnostic>;
}

impl AssemblyExt for Structure {
    fn assembly_set(&self) -> Option<&AssemblySet> {
        self.extensions().get(ASSEMBLIES_EXTENSION)
    }

    fn assembly(&self, id: &str) -> Result<AssemblyView, Diagnostic> {
        let assemblies = self
            .assembly_set()
            .ok_or_else(|| Diagnostic::new(Code::E6002).with_context("assembly", id))?;
        AssemblyView::new(self, assemblies, id)
    }
}

fn resolve_chains(
    structure: &Structure,
    asym_ids: &[Box<str>],
) -> Result<Vec<(ChainIndex, Range<u32>)>, Diagnostic> {
    let mut chains = Vec::with_capacity(asym_ids.len());
    for asym_id in asym_ids {
        let wanted = structure.data().dictionary.get(asym_id).ok_or_else(|| {
            Diagnostic::new(Code::E6013).with_context("label_asym_id", asym_id.as_ref())
        })?;
        let chain = structure
            .data()
            .topology
            .chains
            .iter()
            .find(|chain| structure.data().topology.chains.label_asym_id(*chain) == Some(wanted))
            .ok_or_else(|| {
                Diagnostic::new(Code::E6013).with_context("label_asym_id", asym_id.as_ref())
            })?;
        chains.push((chain, chain_atoms(structure, chain)));
    }
    Ok(chains)
}

fn chain_atoms(structure: &Structure, chain: ChainIndex) -> Range<u32> {
    let Some(residues) = structure.data().topology.chains.residues(chain) else {
        return 0..0;
    };
    let mut atoms = None::<Range<u32>>;
    for residue in residues {
        let Some(range) = structure
            .data()
            .topology
            .residues
            .atoms(molframe_core::ResidueIndex::new(residue))
        else {
            continue;
        };
        match &mut atoms {
            Some(atoms) => atoms.end = range.end,
            None => atoms = Some(range),
        }
    }
    match atoms {
        Some(atoms) => atoms,
        None => 0..0,
    }
}

fn compose(assemblies: &AssemblySet, identifiers: &[&str]) -> Result<Rigid, Diagnostic> {
    let mut transform = Rigid::IDENTITY;
    for identifier in identifiers.iter().rev() {
        let operator = assemblies
            .operator(identifier)
            .ok_or_else(|| Diagnostic::new(Code::E6013).with_context("operator", *identifier))?;
        transform = transform.then(&operator.transform);
    }
    Ok(transform)
}

fn first_unknown_operator(
    assemblies: &AssemblySet,
    expression: &crate::OperExpression,
) -> Diagnostic {
    let identifier = expression
        .factors()
        .iter()
        .flat_map(|factor| factor.iter())
        .find(|identifier| assemblies.operator(identifier).is_none())
        .map_or("(unknown)", Box::as_ref);
    Diagnostic::new(Code::E6013).with_context("operator", identifier)
}

fn instance_limit(assembly: &str, limit: usize) -> Diagnostic {
    Diagnostic::new(Code::E6011)
        .with_context("assembly", assembly)
        .with_context("limit", limit.to_string())
}

#[cfg(test)]
#[path = "view_tests.rs"]
pub(crate) mod tests;
