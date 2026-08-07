//! Chemistry-aware shorthands with deterministic standard-component fallbacks.

use crate::ast::Macro;
use crate::predicate::{AtomContext, scan};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_core::topology::PolymerKind;
use std::collections::BTreeSet;

pub(super) fn macro_selection(
    structure: &Structure,
    universe: &AtomSelection,
    macro_name: Macro,
) -> Result<AtomSelection, Diagnostic> {
    if macro_name == Macro::Aromatic {
        if !structure.data().bonds.is_available() {
            return Err(
                Diagnostic::new(Code::E4003).with_context("required", "aromatic bond annotations")
            );
        }
        let atoms: BTreeSet<_> = structure
            .data()
            .bonds
            .iter()
            .filter(|bond| bond.order == pdbiox_core::BondOrder::Aromatic)
            .flat_map(|bond| [bond.atom_a, bond.atom_b])
            .collect();
        return Ok(scan(structure, universe, |context| {
            atoms.contains(&context.atom.index())
        }));
    }
    Ok(scan(structure, universe, |context| {
        macro_matches(context, macro_name)
    }))
}

fn macro_matches(context: AtomContext<'_>, macro_name: Macro) -> bool {
    let residue_name = context
        .atom
        .component_name()
        .or_else(|| context.residue.name())
        .map_or("", |name| name);
    let atom_name = context.atom.name().map_or("", |name| name);
    let protein =
        context.chain.polymer_kind() == PolymerKind::Protein || AMINO_ACIDS.contains(&residue_name);
    let nucleic = context.chain.polymer_kind().is_nucleic() || NUCLEOTIDES.contains(&residue_name);
    let water = WATER.contains(&residue_name);
    let hydrogen = context
        .atom
        .element()
        .is_some_and(pdbiox_core::element::Element::is_hydrogen);
    let ion = context.residue.is_het() && context.residue.atoms().count() == 1 && !water;
    match macro_name {
        Macro::Protein => protein,
        Macro::Backbone => protein && PROTEIN_BACKBONE.contains(&atom_name),
        Macro::Sidechain => protein && !PROTEIN_BACKBONE.contains(&atom_name) && !hydrogen,
        Macro::Nucleic => nucleic,
        Macro::NucleicBackbone => nucleic && NUCLEIC_BACKBONE.contains(&atom_name),
        Macro::NucleicBase => {
            nucleic && !NUCLEIC_BACKBONE.contains(&atom_name) && !NUCLEIC_SUGAR.contains(&atom_name)
        }
        Macro::NucleicSugar => nucleic && NUCLEIC_SUGAR.contains(&atom_name),
        Macro::Water => water,
        Macro::Ion => ion,
        Macro::Lipid => LIPIDS.contains(&residue_name),
        Macro::Saccharide => {
            context.chain.polymer_kind() == PolymerKind::Saccharide
                || SACCHARIDES.contains(&residue_name)
        }
        Macro::Hetero => context.residue.is_het(),
        Macro::Hydrogen => hydrogen,
        Macro::Heavy => !hydrogen,
        Macro::Polymer => context.chain.polymer_kind().is_polymer() || !context.residue.is_het(),
        Macro::Ligand => context.residue.is_het() && !water && !ion,
        Macro::Aromatic => false,
    }
}

const AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
];
const NUCLEOTIDES: &[&str] = &["A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DI"];
const WATER: &[&str] = &["HOH", "WAT", "H2O", "DOD"];
const PROTEIN_BACKBONE: &[&str] = &["N", "CA", "C", "O"];
const NUCLEIC_BACKBONE: &[&str] = &["P", "OP1", "OP2", "O5'", "C5'", "C4'", "C3'", "O3'"];
const NUCLEIC_SUGAR: &[&str] = &["C1'", "C2'", "C3'", "C4'", "O4'"];
const LIPIDS: &[&str] = &["POPC", "POPE", "DPPC", "DOPC", "CHOL"];
const SACCHARIDES: &[&str] = &["NAG", "BMA", "MAN", "GLC", "GAL", "FUC"];
