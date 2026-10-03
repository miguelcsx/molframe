//! Structure-aligned charges from trusted file columns or a linked CCD graph.
//!
//! Absent CCD hydrogens participate in PEOE, then their charge is folded into
//! their observed parent atom (a united-atom projection). Missing non-leaving
//! heavy atoms are an error. Inter-residue topology joins the calculation;
//! absent leaving atoms attached to a linked atom are removed before typing.

use super::perceive::component_inputs_retained;
use super::{PeoeAtom, PeoeBond, PeoeError, PeoeOptions, peoe_charges};
use crate::{Component, ComponentProvider};
use molframe_core::annotation::{AtomAnnotation, PARTIAL_CHARGE_ANNOTATION};
use molframe_core::contract::DictionaryVersion;
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use molframe_core::{AtomIndex, Diagnostic, Presence};
use std::collections::BTreeMap;

/// Scientific origin of a structure-aligned charge column.
#[derive(Clone, Debug, PartialEq)]
pub enum ChargeSource {
    /// Every atom has a finite, present source-file charge.
    File,
    /// PEOE on the linked component graph, with implicit H charge projected.
    Peoe {
        /// Exact component dictionary used to resolve chemistry.
        dictionary: DictionaryVersion,
        /// Explicit numerical settings.
        options: PeoeOptions,
    },
}

/// One charge in elementary-charge units per structure atom, in atom order.
#[derive(Clone, Debug, PartialEq)]
pub struct PartialCharges {
    /// Charge values aligned with deposited atom indices.
    pub values: Vec<f64>,
    /// Source and calculation policy; file and calculated values are not mixed.
    pub source: ChargeSource,
}

/// Chemistry could not be resolved without inventing charge data.
#[derive(Clone, Debug, PartialEq)]
pub enum PartialChargeError {
    /// Provider storage could not be read.
    Provider(Diagnostic),
    /// Missing component definition or component identifier.
    UnknownComponent(Box<str>),
    /// Observed atom is absent from its component, or has a different element.
    UnknownAtom {
        /// Component identifier.
        component: Box<str>,
        /// Component-local atom identifier.
        atom: Box<str>,
    },
    /// An absent atom cannot be projected onto one observed hydrogen parent.
    MissingAtom {
        /// Component identifier.
        component: Box<str>,
        /// Component-local atom identifier.
        atom: Box<str>,
    },
    /// Alternate conformers share a component atom name and require selection.
    AmbiguousAtom {
        /// Component identifier.
        component: Box<str>,
        /// Repeated component-local atom identifier.
        atom: Box<str>,
    },
    /// Connectivity is unavailable, so polymer linkage cannot be resolved.
    MissingTopology,
    /// A source charge column is malformed or incomplete.
    InvalidFileCharges,
    /// The canonical PEOE kernel rejected component chemistry or options.
    Peoe(PeoeError),
}

impl std::fmt::Display for PartialChargeError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "partial charges cannot resolve chemistry: {self:?}"
        )
    }
}
impl std::error::Error for PartialChargeError {}

/// Calculates charges in stable atom order, preferring a complete finite file
/// column. No source column is mixed with calculated data. A partial/malformed
/// file column is rejected rather than silently replaced.
///
/// The provider makes component and protonation assumptions explicit. File
/// charges do not require dictionary access. PEOE uses all explicit covalent
/// inter-residue bonds and includes CCD hydrogens absent from the coordinates.
///
/// # Errors
///
/// Returns structured errors for unavailable topology, unresolved components,
/// unsupported chemistry, absent heavy atoms, or invalid source charges.
///
/// # Examples
///
/// ```
/// use molframe_chem::{ComponentProvider, PartialCharges, PartialChargeError, PeoeOptions, partial_charges};
/// use molframe_core::structure::Structure;
/// fn charge_column(structure: &Structure, dictionary: &dyn ComponentProvider)
///     -> Result<PartialCharges, PartialChargeError> {
///     partial_charges(structure, dictionary, PeoeOptions::default())
/// }
/// ```
pub fn partial_charges(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    options: PeoeOptions,
) -> Result<PartialCharges, PartialChargeError> {
    if let Some(annotation) = structure.annotations().get(PARTIAL_CHARGE_ANNOTATION) {
        let AtomAnnotation::Real(column) = annotation else {
            return Err(PartialChargeError::InvalidFileCharges);
        };
        if column.len() != structure.atom_count()
            || (0..column.len()).any(|atom| column.presence(atom) != Presence::Present)
            || column.values().iter().any(|charge| !charge.is_finite())
        {
            return Err(PartialChargeError::InvalidFileCharges);
        }
        return Ok(PartialCharges {
            values: column.values().to_vec(),
            source: ChargeSource::File,
        });
    }
    if !structure.data().bonds.is_available() {
        return Err(PartialChargeError::MissingTopology);
    }
    let mut graph = ChargeGraph {
        atoms: Vec::new(),
        bonds: Vec::new(),
        observed: vec![None; structure.atom_count() as usize],
        projected: Vec::new(),
    };
    let mut components = BTreeMap::new();
    for residue in structure.data().residues() {
        let name = match residue.name() {
            Some(name) => name,
            None => "",
        };
        if !components.contains_key(name) {
            let component = provider
                .get(name)
                .map_err(PartialChargeError::Provider)?
                .ok_or_else(|| PartialChargeError::UnknownComponent(name.into()))?;
            components.insert(name, component);
        }
        let component = &components[name];
        graph.add_residue(structure, residue, component)?;
    }
    for bond in structure.data().bonds.iter() {
        let first =
            graph.observed[bond.atom_a.as_usize()].ok_or(PartialChargeError::MissingTopology)?;
        let second =
            graph.observed[bond.atom_b.as_usize()].ok_or(PartialChargeError::MissingTopology)?;
        graph.bonds.push(PeoeBond {
            atom_a: first.min(second),
            atom_b: first.max(second),
        });
    }
    graph
        .bonds
        .sort_unstable_by_key(|bond| (bond.atom_a, bond.atom_b));
    graph.bonds.dedup_by_key(|bond| (bond.atom_a, bond.atom_b));
    let calculated =
        peoe_charges(&graph.atoms, &graph.bonds, options).map_err(PartialChargeError::Peoe)?;
    let mut values = vec![0.0; structure.atom_count() as usize];
    for (charge, parent) in calculated.into_iter().zip(graph.projected) {
        values[parent.as_usize()] += charge;
    }
    Ok(PartialCharges {
        values,
        source: ChargeSource::Peoe {
            dictionary: provider.version().clone(),
            options,
        },
    })
}

struct ChargeGraph {
    atoms: Vec<PeoeAtom>,
    bonds: Vec<PeoeBond>,
    observed: Vec<Option<usize>>,
    projected: Vec<AtomIndex>,
}

impl ChargeGraph {
    fn add_residue(
        &mut self,
        structure: &Structure,
        residue: ResidueRef<'_>,
        component: &Component,
    ) -> Result<(), PartialChargeError> {
        let mut observed = BTreeMap::new();
        for atom in residue.atoms() {
            let name = atom.name().ok_or_else(|| PartialChargeError::UnknownAtom {
                component: component.id.clone(),
                atom: "".into(),
            })?;
            if observed.insert(name, atom).is_some() {
                return Err(PartialChargeError::AmbiguousAtom {
                    component: component.id.clone(),
                    atom: name.into(),
                });
            }
            if component
                .atom(name)
                .is_none_or(|expected| Some(expected.element) != atom.element())
            {
                return Err(PartialChargeError::UnknownAtom {
                    component: component.id.clone(),
                    atom: name.into(),
                });
            }
        }
        let retained = retained_atoms(structure, residue, component, &observed);
        let (atoms, bonds) = component_inputs_retained(component, Some(&retained))
            .map_err(PartialChargeError::Peoe)?;
        let offset = self.atoms.len();
        for (index, atom) in component
            .atoms
            .iter()
            .zip(&retained)
            .filter_map(|(atom, keep)| keep.then_some(atom))
            .enumerate()
        {
            let parent = match observed.get(atom.name.as_ref()) {
                Some(observed) => {
                    self.observed[observed.index().as_usize()] = Some(offset + index);
                    observed.index()
                }
                None if atom.element.is_hydrogen() => {
                    let mut parents = component.bonds.iter().filter_map(|bond| {
                        if bond.atom_a == atom.name {
                            observed.get(bond.atom_b.as_ref())
                        } else if bond.atom_b == atom.name {
                            observed.get(bond.atom_a.as_ref())
                        } else {
                            None
                        }
                    });
                    match (parents.next(), parents.next()) {
                        (Some(parent), None) => parent.index(),
                        _ => {
                            return Err(PartialChargeError::MissingAtom {
                                component: component.id.clone(),
                                atom: atom.name.clone(),
                            });
                        }
                    }
                }
                None => {
                    return Err(PartialChargeError::MissingAtom {
                        component: component.id.clone(),
                        atom: atom.name.clone(),
                    });
                }
            };
            self.projected.push(parent);
        }
        self.atoms.extend(atoms);
        self.bonds.extend(bonds.into_iter().map(|bond| PeoeBond {
            atom_a: offset + bond.atom_a.min(bond.atom_b),
            atom_b: offset + bond.atom_a.max(bond.atom_b),
        }));
        Ok(())
    }
}

fn retained_atoms(
    structure: &Structure,
    residue: ResidueRef<'_>,
    component: &Component,
    observed: &BTreeMap<&str, AtomRef<'_>>,
) -> Vec<bool> {
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    let linked: Vec<_> = observed
        .iter()
        .filter_map(|(name, atom)| {
            adjacency
                .neighbours(atom.index())
                .iter()
                .any(|neighbour| {
                    structure
                        .atom(*neighbour)
                        .and_then(AtomRef::residue)
                        .is_some_and(|other| other.index() != residue.index())
                })
                .then_some(*name)
        })
        .collect();
    let names: BTreeMap<_, _> = component
        .atoms
        .iter()
        .enumerate()
        .map(|(index, atom)| (atom.name.as_ref(), index))
        .collect();
    let mut retained = vec![true; component.atoms.len()];
    loop {
        let mut changed = false;
        for bond in &*component.bonds {
            for (name, other) in [(&bond.atom_a, &bond.atom_b), (&bond.atom_b, &bond.atom_a)] {
                let (Some(&index), Some(&parent)) =
                    (names.get(name.as_ref()), names.get(other.as_ref()))
                else {
                    continue;
                };
                let atom = &component.atoms[index];
                if retained[index]
                    && !observed.contains_key(name.as_ref())
                    && ((atom.leaving && linked.contains(&other.as_ref()))
                        || (!retained[parent] && (atom.leaving || atom.element.is_hydrogen())))
                {
                    retained[index] = false;
                    changed = true;
                }
            }
        }
        if !changed {
            return retained;
        }
    }
}

#[cfg(test)]
#[path = "structure_tests.rs"]
mod tests;
