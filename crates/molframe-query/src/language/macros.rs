//! Chemistry-aware shorthands over explicit structure annotations.

use crate::ast::Macro;
use crate::predicate::{AtomContext, scan};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use std::collections::BTreeSet;

pub(crate) fn macro_selection(
    structure: &Structure,
    universe: &AtomSelection,
    macro_name: Macro,
    _warnings: &mut Vec<Diagnostic>,
) -> Result<AtomSelection, Diagnostic> {
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
    if chemistry_macro(macro_name)
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

const fn chemistry_macro(macro_name: Macro) -> bool {
    matches!(
        macro_name,
        Macro::Protein
            | Macro::Backbone
            | Macro::Sidechain
            | Macro::Nucleic
            | Macro::NucleicBackbone
            | Macro::NucleicBase
            | Macro::NucleicSugar
            | Macro::Water
            | Macro::Ion
            | Macro::Lipid
            | Macro::Saccharide
            | Macro::Ligand
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
    let component_kind = crate::annotation::component_kind(structure, context.atom.index().get());
    let polymer_role = crate::annotation::polymer_atom_role(structure, context.atom.index().get());
    let protein = component_kind == Some(molframe_chem::ComponentKind::AminoAcid);
    let nucleic = component_kind == Some(molframe_chem::ComponentKind::Nucleotide);
    let water = component_kind == Some(molframe_chem::ComponentKind::Solvent);
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
        Macro::Polymer => context.chain.polymer_kind().is_polymer(),
        Macro::Ligand => component_kind == Some(molframe_chem::ComponentKind::NonPolymer),
        Macro::Aromatic => false,
    }
}
