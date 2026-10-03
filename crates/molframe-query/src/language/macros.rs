//! Chemistry-aware shorthands over explicit structure annotations and topology.
//! Hydrogen polarity scans the universe and its hydrogen incident edges, with
//! O(N + E) time and storage on the first use of the shared CSR bond adjacency.

use crate::ast::Macro;
use crate::predicate::{AtomContext, scan};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::element::Element;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{AtomRef, Structure};
use std::collections::BTreeSet;

pub(crate) fn macro_selection(
    structure: &Structure,
    universe: &AtomSelection,
    macro_name: Macro,
    _warnings: &mut Vec<Diagnostic>,
) -> Result<AtomSelection, Diagnostic> {
    if matches!(macro_name, Macro::PolarHydrogen | Macro::NonpolarHydrogen) {
        if !structure.data().bonds.is_available() {
            return Err(Diagnostic::new(Code::E4003).with_context("required", "bond topology"));
        }
        let adjacency = structure.data().bonds.adjacency(structure.atom_count());
        return Ok(scan(structure, universe, |context| {
            if !context.atom.element().is_some_and(Element::is_hydrogen) {
                return false;
            }
            // The universe limits output, not connectivity: a conjunction can
            // exclude the heavy parent before this macro is evaluated.
            let polar = adjacency
                .neighbours(context.atom.index())
                .iter()
                .any(|&neighbour| {
                    matches!(
                        structure.atom(neighbour).and_then(AtomRef::element),
                        Some(Element::NITROGEN | Element::OXYGEN | Element::SULFUR)
                    )
                });
            polar == (macro_name == Macro::PolarHydrogen)
        }));
    }
    if macro_name == Macro::Aromatic {
        if structure
            .annotations()
            .get(molframe_core::AROMATIC_ATOM_ANNOTATION)
            .is_some()
        {
            return Ok(scan(structure, universe, |context| {
                crate::annotation::boolean(
                    structure,
                    context.atom.index().get(),
                    molframe_core::AROMATIC_ATOM_ANNOTATION,
                ) == Some(true)
            }));
        }
        if !structure.data().bonds.is_available() {
            return Err(
                Diagnostic::new(Code::E4003).with_context("required", "aromatic bond annotations")
            );
        }
        let atoms: BTreeSet<_> = structure
            .data()
            .bonds
            .iter()
            .filter(|bond| bond.order == molframe_core::BondOrder::Aromatic)
            .flat_map(|bond| [bond.atom_a, bond.atom_b])
            .collect();
        return Ok(scan(structure, universe, |context| {
            atoms.contains(&context.atom.index())
        }));
    }
    if annotation_only_macro(macro_name)
        && structure
            .annotations()
            .get(molframe_core::COMPONENT_KIND_ANNOTATION)
            .is_none()
    {
        return Err(
            Diagnostic::new(Code::E4003).with_context("required", "CCD component annotations")
        );
    }
    if polymer_role_macro(macro_name)
        && structure
            .annotations()
            .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION)
            .is_none()
    {
        return Err(Diagnostic::new(Code::E4003)
            .with_context("required", "explicit CCD polymer atom-role annotations"));
    }
    Ok(scan(structure, universe, |context| {
        macro_matches(structure, context, macro_name)
    }))
}

const fn polymer_role_macro(macro_name: Macro) -> bool {
    matches!(
        macro_name,
        Macro::Backbone
            | Macro::Sidechain
            | Macro::NucleicBackbone
            | Macro::NucleicBase
            | Macro::NucleicSugar
    )
}

const fn annotation_only_macro(macro_name: Macro) -> bool {
    matches!(
        macro_name,
        Macro::Backbone
            | Macro::Sidechain
            | Macro::NucleicBackbone
            | Macro::NucleicBase
            | Macro::NucleicSugar
            | Macro::Ion
            | Macro::Lipid
            | Macro::Saccharide
    )
}

pub(crate) fn chirality_selection(
    structure: &Structure,
    universe: &AtomSelection,
    configuration: &str,
) -> Result<AtomSelection, Diagnostic> {
    let Some(molframe_core::AtomAnnotation::Symbol(_)) = structure
        .annotations()
        .get(molframe_core::STEREO_CONFIGURATION_ANNOTATION)
    else {
        return Err(Diagnostic::new(Code::E4003).with_context("required", "CCD stereochemistry"));
    };
    Ok(scan(structure, universe, |context| {
        crate::annotation::symbol(
            structure,
            context.atom.index().get(),
            molframe_core::STEREO_CONFIGURATION_ANNOTATION,
        )
        .and_then(|symbol| structure.resolve(symbol))
        .is_some_and(|observed| observed.eq_ignore_ascii_case(configuration))
    }))
}

fn macro_matches(structure: &Structure, context: AtomContext<'_>, macro_name: Macro) -> bool {
    use molframe_core::SecondaryStructure as Ss;
    let secondary = match structure
        .secondary_structure()
        .get(context.residue.index().as_usize())
    {
        Some(state) => *state,
        None => Ss::Unknown,
    };
    let component_kind = crate::annotation::component_kind(structure, context.atom.index().get());
    let polymer_role = crate::annotation::polymer_atom_role(structure, context.atom.index().get());
    let entity_kind = context
        .chain
        .entity()
        .and_then(|entity| structure.data().topology.entities.kind(entity));
    let protein = component_kind == Some(molframe_chem::ComponentKind::AminoAcid)
        || matches!(
            context.chain.polymer_kind(),
            molframe_core::PolymerKind::Protein
        );
    let nucleic = component_kind == Some(molframe_chem::ComponentKind::Nucleotide)
        || context.chain.polymer_kind().is_nucleic();
    let water = component_kind == Some(molframe_chem::ComponentKind::Solvent)
        || entity_kind == Some(molframe_core::EntityKind::Water);
    let hydrogen = context
        .atom
        .element()
        .is_some_and(molframe_core::element::Element::is_hydrogen);
    let ion = component_kind == Some(molframe_chem::ComponentKind::Ion);
    match macro_name {
        Macro::Protein => protein,
        Macro::Backbone => {
            protein
                && polymer_role.is_some_and(|role| {
                    role.intersects(molframe_chem::PolymerAtomRole::PROTEIN_BACKBONE)
                })
        }
        Macro::Sidechain => {
            protein
                && polymer_role.is_some_and(|role| {
                    role.intersects(molframe_chem::PolymerAtomRole::PROTEIN_SIDECHAIN)
                })
                && !hydrogen
        }
        Macro::Nucleic => nucleic,
        Macro::NucleicBackbone => {
            nucleic
                && polymer_role.is_some_and(|role| {
                    role.intersects(molframe_chem::PolymerAtomRole::NUCLEIC_BACKBONE)
                })
        }
        Macro::NucleicBase => {
            nucleic
                && polymer_role.is_some_and(|role| {
                    role.intersects(molframe_chem::PolymerAtomRole::NUCLEIC_BASE_GROUP)
                })
        }
        Macro::NucleicSugar => {
            nucleic
                && polymer_role.is_some_and(|role| {
                    role.intersects(molframe_chem::PolymerAtomRole::NUCLEIC_SUGAR)
                })
        }
        Macro::Water => water,
        Macro::Ion => ion,
        Macro::Lipid => component_kind == Some(molframe_chem::ComponentKind::Lipid),
        Macro::Saccharide => component_kind == Some(molframe_chem::ComponentKind::Saccharide),
        Macro::Hetero => context.residue.is_het(),
        Macro::Hydrogen => hydrogen,
        Macro::Heavy => !hydrogen,
        Macro::Polymer => context.chain.polymer_kind().is_polymer() || protein || nucleic,
        Macro::Ligand => {
            component_kind == Some(molframe_chem::ComponentKind::NonPolymer)
                || entity_kind == Some(molframe_core::EntityKind::NonPolymer)
        }
        Macro::Helix => secondary.is_helix(),
        Macro::Strand => secondary.is_strand(),
        Macro::Sheet => secondary.is_sheet_like(),
        Macro::AlphaHelix => secondary == Ss::AlphaHelix,
        Macro::Helix310 => secondary == Ss::ThreeTenHelix,
        Macro::PiHelix => secondary == Ss::PiHelix,
        Macro::Polyproline => secondary == Ss::PolyProline,
        Macro::Bridge => secondary == Ss::BetaBridge,
        Macro::Turn => secondary == Ss::Turn,
        Macro::Bend => secondary == Ss::Bend,
        Macro::Coil => secondary == Ss::Coil,
        Macro::Aromatic | Macro::PolarHydrogen | Macro::NonpolarHydrogen => false,
    }
}

#[cfg(test)]
#[path = "macros_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod secondary_tests;
