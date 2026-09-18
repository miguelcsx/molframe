//! Pluggable deterministic component providers and CCD lowering.

mod validation;

use crate::model::{Component, ComponentAtom, ComponentBond, ComponentKind, StereoConfiguration};
use molframe_cif::{CifValue, DataBlock, Document};
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::InputBuffer;
use molframe_core::{BondOrder, Element};
use num_traits::ToPrimitive;
use std::collections::BTreeMap;
use std::sync::Arc;
use validation::validate_component_topology;

type ComponentCoordinates = Arc<[[f32; 3]]>;
type LoweringResult<T> = Result<T, Vec<Diagnostic>>;

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
    /// Builds an index without replacing duplicate component definitions.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when two supplied components share an identifier.
    pub fn new(
        version: DictionaryVersion,
        components: impl IntoIterator<Item = Component>,
    ) -> Result<Self, Diagnostic> {
        let mut indexed = BTreeMap::new();
        for component in components {
            let id = component.id.clone();
            if indexed.insert(id.clone(), Arc::new(component)).is_some() {
                return Err(Diagnostic::new(Code::E2005)
                    .in_category("chem_comp")
                    .with_context("id", id.to_string()));
            }
        }
        Ok(Self {
            version,
            components: indexed,
        })
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
    /// Returns every missing, malformed, or duplicate component field without
    /// constructing a partial dictionary.
    pub fn from_document(
        document: &Document,
        version: DictionaryVersion,
    ) -> Result<Self, Vec<Diagnostic>> {
        let mut components = Vec::new();
        let mut findings = Vec::new();
        let mut identifiers = std::collections::BTreeSet::new();
        for block in document.blocks() {
            match component(block) {
                Ok(component) if identifiers.insert(component.id.clone()) => {
                    components.push(component);
                }
                Ok(component) => findings.push(
                    Diagnostic::new(Code::E2005)
                        .in_category("chem_comp")
                        .with_context("id", component.id.to_string())
                        .with_context("block", block.name()),
                ),
                Err(block_findings) => findings.extend(block_findings),
            }
        }
        if findings.is_empty() {
            MemoryProvider::new(version, components)
                .map(Self)
                .map_err(one_finding)
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
    let (document, findings) = molframe_cif::parse(input)?;
    CifProvider::from_document(&document, version).map(|provider| (provider, findings))
}

fn component(block: &DataBlock) -> LoweringResult<Component> {
    let Some(metadata) = block.category("chem_comp") else {
        return Err(vec![
            Diagnostic::new(Code::E2001)
                .with_context("block", block.name())
                .in_category("chem_comp"),
        ]);
    };
    let id = required_identifier(metadata, "id", 0).map_err(one_finding)?;
    let name = required_identifier(metadata, "name", 0).map_err(one_finding)?;
    let component_type = required_text(metadata, "type", 0).map_err(one_finding)?;
    let atoms = atoms(block)?;
    let bonds = bonds(block)?;
    validate_component_topology(&atoms, &bonds)?;
    let component_kind = classify_component(&component_type, &atoms, &bonds);
    let ideal_coordinates = coordinates(block, true, atoms.len())?;
    let model_coordinates = coordinates(block, false, atoms.len())?;
    Ok(Component {
        id,
        name,
        kind: component_kind,
        parent: metadata
            .identifier("mon_nstd_parent_comp_id", 0)
            .map(|value| value.into_owned().into()),
        one_letter_code: metadata
            .identifier("one_letter_code", 0)
            .and_then(|value| one_letter_code(value.as_ref())),
        formula: metadata
            .identifier("formula", 0)
            .map(|value| value.into_owned().into()),
        ideal_coordinates,
        model_coordinates,
        atoms: atoms.into(),
        bonds: bonds.into(),
    })
}

fn classify_component(
    component_type: &str,
    atoms: &[ComponentAtom],
    bonds: &[ComponentBond],
) -> ComponentKind {
    let broad = kind(Some(component_type));
    if broad == ComponentKind::NonPolymer
        && atoms.len() == 1
        && bonds.is_empty()
        && atoms[0].charge != 0
    {
        ComponentKind::Ion
    } else {
        broad
    }
}

fn one_letter_code(value: &str) -> Option<u8> {
    let bytes = value.as_bytes();
    (bytes.len() == 1 && bytes[0].is_ascii()).then_some(bytes[0])
}

fn atoms(block: &DataBlock) -> LoweringResult<Vec<ComponentAtom>> {
    let Some(category) = block.category("chem_comp_atom") else {
        return Err(vec![
            Diagnostic::new(Code::E2001)
                .with_context("block", block.name())
                .in_category("chem_comp_atom"),
        ]);
    };
    if category.row_count() == 0 {
        return Err(vec![
            Diagnostic::new(Code::E2002)
                .with_context("item", "atom_id")
                .with_context("row", "0")
                .in_category(category.name()),
        ]);
    }
    collect_rows(category, |row| {
        let symbol = required_text(category, "type_symbol", row)?;
        let element = Element::from_symbol(&symbol)
            .ok_or_else(|| invalid_value(category, "type_symbol", row))?;
        let charge_value = required_integer(category, "charge", row)?;
        let charge =
            i8::try_from(charge_value).map_err(|_| invalid_value(category, "charge", row))?;
        Ok(ComponentAtom {
            name: required_identifier(category, "atom_id", row)?,
            alternate_name: category
                .identifier("alt_atom_id", row)
                .map(|value| value.into_owned().into()),
            element,
            charge,
            aromatic: required_flag(category, "pdbx_aromatic_flag", row)?,
            leaving: required_flag(category, "pdbx_leaving_atom_flag", row)?,
            stereo: optional_stereo(category, row)?,
        })
    })
}

fn bonds(block: &DataBlock) -> LoweringResult<Vec<ComponentBond>> {
    let Some(category) = block.category("chem_comp_bond") else {
        return Ok(Vec::new());
    };
    collect_rows(category, |row| {
        let aromatic = required_flag(category, "pdbx_aromatic_flag", row)?;
        Ok(ComponentBond {
            atom_a: required_identifier(category, "atom_id_1", row)?,
            atom_b: required_identifier(category, "atom_id_2", row)?,
            order: bond_order(category, row, aromatic)?,
            aromatic,
            stereo: optional_stereo(category, row)?,
        })
    })
}

fn coordinates(
    block: &DataBlock,
    ideal: bool,
    expected: usize,
) -> LoweringResult<Option<ComponentCoordinates>> {
    let Some(category) = block.category("chem_comp_atom") else {
        return Ok(None);
    };
    let items = if ideal {
        [
            "pdbx_model_Cartn_x_ideal",
            "pdbx_model_Cartn_y_ideal",
            "pdbx_model_Cartn_z_ideal",
        ]
    } else {
        ["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"]
    };
    if items.iter().all(|item| category.column(item).is_none()) {
        return Ok(None);
    }
    let mut positions = Vec::with_capacity(expected);
    for row in 0..category.row_count() {
        let mut position = [0.0; 3];
        for (axis, item) in items.iter().enumerate() {
            let value = required_number(category, item, row).map_err(one_finding)?;
            if value < f64::from(f32::MIN) || value > f64::from(f32::MAX) {
                return Err(vec![invalid_value(category, item, row)]);
            }
            let Some(value) = value.to_f32() else {
                return Err(vec![invalid_value(category, item, row)]);
            };
            position[axis] = value;
        }
        positions.push(position);
    }
    if positions.len() == expected {
        Ok(Some(positions.into()))
    } else {
        Err(vec![
            Diagnostic::new(Code::E3001)
                .in_category(category.name())
                .with_context("expected", expected.to_string())
                .with_context("observed", positions.len().to_string()),
        ])
    }
}

fn kind(value: Option<&str>) -> ComponentKind {
    const PEPTIDE: &str = "PEPTIDE";
    const DNA: &str = "DNA";
    const RNA: &str = "RNA";
    const SACCHARIDE: &str = "SACCHARIDE";
    const LIPID: &str = "LIPID";
    const WATER: &str = "WATER";
    const NON_POLYMER: &str = "NON-POLYMER";

    let Some(text) = value else {
        return ComponentKind::Unknown;
    };
    let normalized = text.trim().to_ascii_uppercase();
    let has_word = |word: &str| {
        normalized
            .split(|character: char| !character.is_ascii_alphanumeric())
            .any(|token| token == word)
    };
    if has_word(PEPTIDE) {
        ComponentKind::AminoAcid
    } else if has_word(DNA) || has_word(RNA) {
        ComponentKind::Nucleotide
    } else if has_word(SACCHARIDE) {
        ComponentKind::Saccharide
    } else if has_word(LIPID) {
        ComponentKind::Lipid
    } else if normalized == WATER {
        ComponentKind::Solvent
    } else if normalized == NON_POLYMER {
        ComponentKind::NonPolymer
    } else {
        ComponentKind::Unknown
    }
}

fn bond_order(
    category: &molframe_cif::Category,
    row: usize,
    aromatic: bool,
) -> Result<BondOrder, Diagnostic> {
    if aromatic {
        return Ok(BondOrder::Aromatic);
    }
    match required_text(category, "value_order", row)?
        .to_ascii_uppercase()
        .as_str()
    {
        "SING" => Ok(BondOrder::Single),
        "DOUB" => Ok(BondOrder::Double),
        "TRIP" => Ok(BondOrder::Triple),
        "QUAD" => Ok(BondOrder::Quadruple),
        _ => Err(invalid_value(category, "value_order", row)),
    }
}

fn optional_stereo(
    category: &molframe_cif::Category,
    row: usize,
) -> Result<Option<StereoConfiguration>, Diagnostic> {
    let item = "pdbx_stereo_config";
    match category.text(item, row) {
        None => Ok(None),
        Some(value) if value.eq_ignore_ascii_case("N") => Ok(None),
        Some(value) if value.eq_ignore_ascii_case("R") => Ok(Some(StereoConfiguration::R)),
        Some(value) if value.eq_ignore_ascii_case("S") => Ok(Some(StereoConfiguration::S)),
        Some(value) if value.eq_ignore_ascii_case("MIXED") => Ok(Some(StereoConfiguration::Mixed)),
        Some(_) => Err(invalid_value(category, item, row)),
    }
}

fn collect_rows<T>(
    category: &molframe_cif::Category,
    parse: impl Fn(usize) -> Result<T, Diagnostic>,
) -> LoweringResult<Vec<T>> {
    let mut values = Vec::with_capacity(category.row_count());
    let mut findings = Vec::new();
    for row in 0..category.row_count() {
        match parse(row) {
            Ok(value) => values.push(value),
            Err(finding) => findings.push(finding),
        }
    }
    if findings.is_empty() {
        Ok(values)
    } else {
        Err(findings)
    }
}

fn required_identifier(
    category: &molframe_cif::Category,
    item: &str,
    row: usize,
) -> Result<Box<str>, Diagnostic> {
    category
        .identifier(item, row)
        .map(|value| value.into_owned().into())
        .ok_or_else(|| missing_item(category, item, row))
}

fn required_text(
    category: &molframe_cif::Category,
    item: &str,
    row: usize,
) -> Result<Box<str>, Diagnostic> {
    category
        .text(item, row)
        .map(Into::into)
        .ok_or_else(|| missing_item(category, item, row))
}

fn required_integer(
    category: &molframe_cif::Category,
    item: &str,
    row: usize,
) -> Result<i64, Diagnostic> {
    let value = category
        .value(item, row)
        .ok_or_else(|| missing_item(category, item, row))?;
    if matches!(value, CifValue::Unknown | CifValue::Inapplicable) {
        return Err(missing_item(category, item, row));
    }
    value
        .as_integer()
        .ok_or_else(|| wrong_type(category, item, row))
}

fn required_number(
    category: &molframe_cif::Category,
    item: &str,
    row: usize,
) -> Result<f64, Diagnostic> {
    let value = category
        .value(item, row)
        .ok_or_else(|| missing_item(category, item, row))?;
    if matches!(value, CifValue::Unknown | CifValue::Inapplicable) {
        return Err(missing_item(category, item, row));
    }
    value
        .as_float()
        .filter(|number| number.is_finite())
        .ok_or_else(|| wrong_type(category, item, row))
}

fn one_finding(finding: Diagnostic) -> Vec<Diagnostic> {
    vec![finding]
}

fn required_flag(
    category: &molframe_cif::Category,
    item: &str,
    row: usize,
) -> Result<bool, Diagnostic> {
    match category.text(item, row) {
        Some(value) if value.eq_ignore_ascii_case("Y") => Ok(true),
        Some(value) if value.eq_ignore_ascii_case("N") => Ok(false),
        Some(_) => Err(invalid_value(category, item, row)),
        None => Err(missing_item(category, item, row)),
    }
}

fn missing_item(category: &molframe_cif::Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2002)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

fn wrong_type(category: &molframe_cif::Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2004)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

fn invalid_value(category: &molframe_cif::Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2003)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

#[cfg(test)]
#[path = "provider_tests.rs"]
mod tests;
