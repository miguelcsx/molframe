//! Pluggable deterministic component providers and CCD lowering.

use crate::model::{Component, ComponentAtom, ComponentBond, ComponentKind, StereoConfiguration};
use pdbiox_cif::{CifValue, DataBlock, Document};
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::InputBuffer;
use pdbiox_core::{BondOrder, Element};
use std::collections::BTreeMap;
use std::sync::Arc;

/// A source of chemical component definitions.
pub trait ComponentProvider: Send + Sync {
    /// Returns a component, or `None` when this provider does not know it.
    ///
    /// # Errors
    ///
    /// Returns a registered finding when provider storage cannot be read.
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic>;

    /// Exact dictionary version supplied by this provider.
    fn version(&self) -> &DictionaryVersion;
}

/// An in-memory deterministic provider.
#[derive(Clone, Debug)]
pub struct MemoryProvider {
    version: DictionaryVersion,
    components: BTreeMap<Box<str>, Arc<Component>>,
}

impl MemoryProvider {
    /// Builds an index. Later duplicate identifiers replace earlier ones.
    #[must_use]
    pub fn new(
        version: DictionaryVersion,
        components: impl IntoIterator<Item = Component>,
    ) -> Self {
        let components = components
            .into_iter()
            .map(|component| (component.id.clone(), Arc::new(component)))
            .collect();
        Self {
            version,
            components,
        }
    }
}

impl ComponentProvider for MemoryProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        Ok(self.components.get(component_id).cloned())
    }

    fn version(&self) -> &DictionaryVersion {
        &self.version
    }
}

/// Component provider lowered from a CCD mmCIF document.
#[derive(Clone, Debug)]
pub struct CifProvider(MemoryProvider);

impl CifProvider {
    /// Lowers every component block in a parsed CCD document.
    ///
    /// # Errors
    ///
    /// Returns one diagnostic per block missing its component definition.
    pub fn from_document(
        document: &Document,
        version: DictionaryVersion,
    ) -> Result<Self, Vec<Diagnostic>> {
        let mut components = Vec::new();
        let mut findings = Vec::new();
        for block in document.blocks() {
            match component(block) {
                Some(component) => components.push(component),
                None => findings.push(
                    Diagnostic::new(Code::E2001)
                        .with_context("block", block.name())
                        .in_category("chem_comp"),
                ),
            }
        }
        if findings.is_empty() {
            Ok(Self(MemoryProvider::new(version, components)))
        } else {
            Err(findings)
        }
    }
}

impl ComponentProvider for CifProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        self.0.get(component_id)
    }

    fn version(&self) -> &DictionaryVersion {
        self.0.version()
    }
}

/// Reads a CCD document and lowers every component through the shared CIF parser.
///
/// # Errors
///
/// Returns parser diagnostics or component-definition errors without producing
/// a partial provider.
pub fn read_ccd(
    input: &InputBuffer,
    version: DictionaryVersion,
) -> Result<(CifProvider, Vec<Diagnostic>), Vec<Diagnostic>> {
    let (document, findings) = pdbiox_cif::parse(input)?;
    CifProvider::from_document(&document, version).map(|provider| (provider, findings))
}

fn component(block: &DataBlock) -> Option<Component> {
    let metadata = block.category("chem_comp")?;
    let id = metadata.identifier("id", 0)?.into_owned().into();
    let name = match metadata.identifier("name", 0) {
        Some(name) => name.into_owned().into(),
        None => block.name().into(),
    };
    let atoms = atoms(block);
    Some(Component {
        id,
        name,
        kind: kind(metadata.text("type", 0)),
        parent: metadata
            .identifier("mon_nstd_parent_comp_id", 0)
            .map(|value| value.into_owned().into()),
        formula: metadata
            .identifier("formula", 0)
            .map(|value| value.into_owned().into()),
        ideal_coordinates: coordinates(block, true, atoms.len()),
        model_coordinates: coordinates(block, false, atoms.len()),
        atoms: atoms.into(),
        bonds: bonds(block).into(),
    })
}

fn atoms(block: &DataBlock) -> Vec<ComponentAtom> {
    let Some(category) = block.category("chem_comp_atom") else {
        return Vec::new();
    };
    (0..category.row_count())
        .filter_map(|row| {
            let element = category
                .text("type_symbol", row)
                .and_then(Element::from_symbol);
            let charge = category
                .value("charge", row)
                .and_then(CifValue::as_integer)
                .and_then(|value| i8::try_from(value).ok());
            Some(ComponentAtom {
                name: category.identifier("atom_id", row)?.into_owned().into(),
                alternate_name: category
                    .identifier("alt_atom_id", row)
                    .map(|value| value.into_owned().into()),
                element: match element {
                    Some(element) => element,
                    None => Element::UNKNOWN,
                },
                charge: match charge {
                    Some(charge) => charge,
                    None => 0,
                },
                aromatic: yes(category.text("pdbx_aromatic_flag", row)),
                leaving: yes(category.text("pdbx_leaving_atom_flag", row)),
                stereo: stereo(category.text("pdbx_stereo_config", row)),
            })
        })
        .collect()
}

fn bonds(block: &DataBlock) -> Vec<ComponentBond> {
    let Some(category) = block.category("chem_comp_bond") else {
        return Vec::new();
    };
    (0..category.row_count())
        .filter_map(|row| {
            let aromatic = yes(category.text("pdbx_aromatic_flag", row));
            Some(ComponentBond {
                atom_a: category.identifier("atom_id_1", row)?.into_owned().into(),
                atom_b: category.identifier("atom_id_2", row)?.into_owned().into(),
                order: bond_order(category.text("value_order", row), aromatic),
                aromatic,
                stereo: stereo(category.text("pdbx_stereo_config", row)),
            })
        })
        .collect()
}

fn coordinates(block: &DataBlock, ideal: bool, expected: usize) -> Option<Arc<[[f32; 3]]>> {
    let category = block.category("chem_comp_atom")?;
    let mut positions = Vec::with_capacity(expected);
    for row in 0..category.row_count() {
        let items = if ideal {
            [
                "pdbx_model_Cartn_x_ideal",
                "pdbx_model_Cartn_y_ideal",
                "pdbx_model_Cartn_z_ideal",
            ]
        } else {
            ["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"]
        };
        let read = |item: &str| {
            category
                .value(item, row)
                .and_then(CifValue::as_float)
                .map(|value| value as f32)
        };
        positions.push([read(items[0])?, read(items[1])?, read(items[2])?]);
    }
    (positions.len() == expected).then(|| positions.into())
}

fn kind(value: Option<&str>) -> ComponentKind {
    let lower = value.map(str::to_ascii_lowercase);
    let Some(text) = lower.as_deref() else {
        return ComponentKind::Unknown;
    };
    if text.contains("peptide") {
        ComponentKind::AminoAcid
    } else if text.contains("dna") || text.contains("rna") {
        ComponentKind::Nucleotide
    } else if text.contains("saccharide") {
        ComponentKind::Saccharide
    } else if text.contains("water") {
        ComponentKind::Solvent
    } else if text.contains("non-polymer") {
        ComponentKind::NonPolymer
    } else {
        ComponentKind::Unknown
    }
}

fn bond_order(value: Option<&str>, aromatic: bool) -> BondOrder {
    if aromatic {
        return BondOrder::Aromatic;
    }
    match value.map(str::to_ascii_uppercase).as_deref() {
        Some("SING") => BondOrder::Single,
        Some("DOUB") => BondOrder::Double,
        Some("TRIP") => BondOrder::Triple,
        Some("QUAD") => BondOrder::Quadruple,
        _ => BondOrder::Unknown,
    }
}

fn stereo(value: Option<&str>) -> Option<StereoConfiguration> {
    match value.map(str::to_ascii_uppercase).as_deref() {
        Some("R") => Some(StereoConfiguration::R),
        Some("S") => Some(StereoConfiguration::S),
        Some("MIXED") => Some(StereoConfiguration::Mixed),
        _ => None,
    }
}

fn yes(value: Option<&str>) -> bool {
    value.is_some_and(|value| value.eq_ignore_ascii_case("Y"))
}

#[cfg(test)]
#[path = "provider_tests.rs"]
mod tests;
