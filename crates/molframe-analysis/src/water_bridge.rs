//! Water-mediated networks assembled from fully perceived hydrogen bonds.

use std::collections::BTreeMap;

use molframe_core::ExecutionContext;
use molframe_core::index::AtomIndex;
use molframe_core::structure::Structure;

use crate::hbond::{
    HydrogenBondError, HydrogenBondOptions, HydrogenBondPolicy, altlocs_compatible,
    heavy_atom_pairs, hydrogen_bonds_with_policy,
};

/// Explicit water-bridge detection policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterBridgeOptions {
    /// Geometry and periodic policy for each constituent hydrogen bond.
    pub hydrogen_bonds: HydrogenBondOptions,
}

/// A water heavy atom bridging two non-water partners.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WaterBridge {
    /// Water donor or acceptor atom.
    pub water: AtomIndex,
    /// Lower-indexed non-water partner.
    pub first: AtomIndex,
    /// Higher-indexed non-water partner.
    pub second: AtomIndex,
}

define_soa_table! {
    /// Native columnar storage for water-mediated interactions.
    pub struct WaterBridgeTable for WaterBridge {
        /// Bridging solvent atom indices.
        water: AtomIndex,
        /// Lower partner atom indices.
        first: AtomIndex,
        /// Higher partner atom indices.
        second: AtomIndex,
    }
}

/// How hydrogen-bond evidence for each water contact is gathered.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WaterBridgeMode {
    /// Every contact must pass the explicit-hydrogen geometry test.
    #[default]
    Strict,
    /// Contacts are donor/acceptor heavy-atom distances only, for structures
    /// whose waters carry no hydrogens. Lower confidence: no angle is checked,
    /// and a pair of acceptors can look like a bridge.
    HeavyAtomOnly,
}

/// Evidence level behind the bridges in a [`WaterBridgeReport`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterBridgeConfidence {
    /// Each constituent contact has an explicit hydrogen and an angle.
    Hydrogen,
    /// Contacts rest on heavy-atom distances alone.
    HeavyAtomOnly,
}

/// Mode and hydrogen-bond refinements for [`water_bridges_with_policy`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WaterBridgePolicy {
    /// Evidence mode.
    pub mode: WaterBridgeMode,
    /// Refinements applied to each constituent hydrogen bond.
    pub hydrogen_bonds: HydrogenBondPolicy,
}

/// Bridges plus how trustworthy and how complete the search was.
#[derive(Clone, Debug)]
pub struct WaterBridgeReport {
    /// The detected bridges.
    pub bridges: WaterBridgeTable,
    /// Evidence level of every bridge in the table.
    pub confidence: WaterBridgeConfidence,
    /// Water molecules that carry no hydrogen atoms at all.
    ///
    /// In strict mode such a water can only be an acceptor, so bridges through
    /// it may be missed; the count says how many waters were not fully
    /// evaluable.
    pub waters_without_hydrogens: usize,
}

const WATER_COMPONENT_IDS: [&str; 7] = ["HOH", "WAT", "DOD", "H2O", "TIP3", "TIP", "SOL"];

/// Whether the atom belongs to a water molecule.
///
/// Water is recognised by its component identifier or by composition (one
/// oxygen and at most two hydrogens, for components not classified as a
/// non-solvent kind). A solvent such as glycerol or DMSO is a
/// solvent but not water, so it never bridges.
fn is_water(structure: &Structure, atom: u32) -> bool {
    let Some(atom) = structure.data().atom(AtomIndex::new(atom)) else {
        return false;
    };
    if atom.component_name().is_some_and(|name| {
        WATER_COMPONENT_IDS
            .iter()
            .any(|id| id.eq_ignore_ascii_case(name.trim()))
    }) {
        return true;
    }
    // Composition only decides for components not already classified as
    // something else: an amino-acid or ligand residue with one oxygen is not
    // water.
    let classified_other = crate::chemistry::component_kind(structure, atom.index().get())
        .is_some_and(|kind| kind != molframe_chem::ComponentKind::Solvent);
    !classified_other && atom.residue().is_some_and(water_composition)
}

fn water_composition(residue: molframe_core::structure::ResidueRef<'_>) -> bool {
    let (mut oxygens, mut hydrogens, mut others) = (0, 0, 0);
    for atom in residue.atoms() {
        match atom.element() {
            Some(element) if element.is_hydrogen() => hydrogens += 1,
            Some(element) if element.atomic_number() == 8 => oxygens += 1,
            _ => others += 1,
        }
    }
    oxygens == 1 && others == 0 && hydrogens <= 2
}

fn waters_without_hydrogens(structure: &Structure) -> usize {
    structure
        .data()
        .residues()
        .filter(|residue| {
            let has_oxygen = residue
                .atoms()
                .any(|atom| atom.element().is_some_and(|e| e.atomic_number() == 8));
            has_oxygen
                && residue
                    .atoms()
                    .next()
                    .is_some_and(|atom| is_water(structure, atom.index().get()))
                && !residue.atoms().any(|atom| {
                    atom.element()
                        .is_some_and(molframe_core::Element::is_hydrogen)
                })
        })
        .count()
}

/// Builds water bridges from two or more oriented hydrogen bonds to one water.
///
/// Component identity and donor/acceptor roles come exclusively from CCD
/// annotations. Only water bridges: other solvents (glycerol, DMSO) are not
/// water and never mediate. Two partners in incompatible alternate locations
/// are never joined. This is the strict, explicit-hydrogen search; see
/// [`water_bridges_with_policy`] for the heavy-atom mode and the report.
///
/// # Errors
///
/// Propagates hydrogen-bond chemistry, geometry, periodic and spatial errors.
pub fn water_bridges(
    structure: &Structure,
    options: WaterBridgeOptions,
    context: &ExecutionContext,
) -> Result<WaterBridgeTable, HydrogenBondError> {
    water_bridges_with_policy(structure, options, WaterBridgePolicy::default(), context)
        .map(|report| report.bridges)
}

/// Builds water bridges under an explicit [`WaterBridgePolicy`].
///
/// # Errors
///
/// Propagates hydrogen-bond chemistry, geometry, periodic and spatial errors.
pub fn water_bridges_with_policy(
    structure: &Structure,
    options: WaterBridgeOptions,
    policy: WaterBridgePolicy,
    context: &ExecutionContext,
) -> Result<WaterBridgeReport, HydrogenBondError> {
    let oriented: Vec<(u32, u32)> = match policy.mode {
        WaterBridgeMode::Strict => hydrogen_bonds_with_policy(
            structure,
            options.hydrogen_bonds,
            policy.hydrogen_bonds,
            context,
        )?
        .iter()
        .map(|bond| (bond.donor.get(), bond.acceptor.get()))
        .collect(),
        WaterBridgeMode::HeavyAtomOnly => heavy_atom_pairs(
            structure,
            options.hydrogen_bonds,
            policy.hydrogen_bonds,
            context,
        )?,
    };
    let mut partners: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (donor, acceptor) in oriented {
        match (is_water(structure, donor), is_water(structure, acceptor)) {
            (true, false) => partners.entry(donor).or_default().push(acceptor),
            (false, true) => partners.entry(acceptor).or_default().push(donor),
            (true, true) | (false, false) => {}
        }
    }
    let mut bridges = Vec::new();
    for (water, mut atoms) in partners {
        atoms.sort_unstable();
        atoms.dedup();
        for first in 0..atoms.len() {
            for second in first + 1..atoms.len() {
                if altlocs_compatible(structure, atoms[first], atoms[second]) {
                    bridges.push(WaterBridge {
                        water: AtomIndex::new(water),
                        first: AtomIndex::new(atoms[first]),
                        second: AtomIndex::new(atoms[second]),
                    });
                }
            }
        }
    }
    Ok(WaterBridgeReport {
        bridges: bridges.into_iter().collect(),
        confidence: match policy.mode {
            WaterBridgeMode::Strict => WaterBridgeConfidence::Hydrogen,
            WaterBridgeMode::HeavyAtomOnly => WaterBridgeConfidence::HeavyAtomOnly,
        },
        waters_without_hydrogens: waters_without_hydrogens(structure),
    })
}

#[cfg(test)]
#[path = "water_bridge_tests.rs"]
mod tests;
