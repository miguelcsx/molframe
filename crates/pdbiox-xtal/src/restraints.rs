//! CCP4/Refmac-style monomer restraint libraries.

use pdbiox_cif::{CifValue, DataBlock, Document};
use pdbiox_core::io::InputBuffer;
use std::collections::{BTreeMap, BTreeSet};

/// Ideal bond-length restraint.
#[derive(Clone, Debug, PartialEq)]
pub struct BondRestraint {
    /// Component-local atom identifiers.
    pub atoms: [Box<str>; 2],
    /// Target distance in ångströms.
    pub target: f64,
    /// Standard uncertainty in ångströms.
    pub sigma: f64,
    /// Deposited chemical restraint type.
    pub kind: Option<Box<str>>,
}

/// Ideal valence-angle restraint.
#[derive(Clone, Debug, PartialEq)]
pub struct AngleRestraint {
    /// Atom identifiers; the middle atom is the vertex.
    pub atoms: [Box<str>; 3],
    /// Target angle in degrees.
    pub target: f64,
    /// Standard uncertainty in degrees.
    pub sigma: f64,
}

/// Ideal torsion restraint.
#[derive(Clone, Debug, PartialEq)]
pub struct TorsionRestraint {
    /// Stable torsion identifier or label.
    pub id: Option<Box<str>>,
    /// Four ordered atom identifiers.
    pub atoms: [Box<str>; 4],
    /// Target dihedral angle in degrees.
    pub target: f64,
    /// Standard uncertainty in degrees.
    pub sigma: f64,
    /// Number of equivalent minima over 360 degrees.
    pub period: u32,
}

/// One atom's membership in a planar group.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneAtomRestraint {
    /// Component-local atom identifier.
    pub atom: Box<str>,
    /// Standard uncertainty of distance from the plane in ångströms.
    pub sigma: f64,
}

/// Planarity restraint grouped by plane identifier.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneRestraint {
    /// Plane identifier.
    pub id: Box<str>,
    /// Atoms expected in the plane.
    pub atoms: Vec<PlaneAtomRestraint>,
}

/// Expected sign of a tetrahedral chiral volume.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChiralVolumeSign {
    /// Positive signed volume.
    Positive,
    /// Negative signed volume.
    Negative,
    /// Either sign is permitted.
    Both,
}

/// Chiral-volume restraint.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChiralRestraint {
    /// Stable chirality identifier.
    pub id: Option<Box<str>>,
    /// Central atom followed by three ordered substituents.
    pub atoms: [Box<str>; 4],
    /// Required volume sign.
    pub sign: ChiralVolumeSign,
}

/// Complete deterministic restraint set for one monomer.
#[derive(Clone, Debug, PartialEq)]
pub struct MonomerRestraints {
    /// Component identifier.
    pub id: Box<str>,
    /// Declared atom identifiers.
    pub atoms: Vec<Box<str>>,
    /// Bond restraints.
    pub bonds: Vec<BondRestraint>,
    /// Angle restraints.
    pub angles: Vec<AngleRestraint>,
    /// Torsion restraints.
    pub torsions: Vec<TorsionRestraint>,
    /// Planar groups.
    pub planes: Vec<PlaneRestraint>,
    /// Chiral-volume restraints.
    pub chirals: Vec<ChiralRestraint>,
}

/// Indexed monomer restraint library.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MonomerLibrary {
    components: BTreeMap<Box<str>, MonomerRestraints>,
}

impl MonomerLibrary {
    /// Retrieves one component case-sensitively.
    #[must_use]
    pub fn get(&self, id: &str) -> Option<&MonomerRestraints> {
        self.components.get(id)
    }

    /// Number of indexed components.
    #[must_use]
    pub fn len(&self) -> usize {
        self.components.len()
    }

    /// Returns true when no component restraint sets were loaded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }

    /// Iterates in stable component-identifier order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &MonomerRestraints)> {
        self.components
            .iter()
            .map(|(id, restraints)| (id.as_ref(), restraints))
    }
}

/// Invalid monomer restraint library.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum RestraintError {
    /// No component restraint blocks were present.
    #[error("monomer library has no component restraint blocks")]
    MissingComponents,
    /// Component identifiers are absent or duplicated.
    #[error("monomer library component identifier is absent or duplicated")]
    ComponentId,
    /// A restraint names an atom absent from its component.
    #[error("monomer restraint references an absent component atom")]
    AtomReference,
    /// A required target, sigma, period, or chirality sign is malformed.
    #[error("monomer restraint contains an invalid target or type")]
    InvalidValue,
}

/// Lowers all component blocks of a CCP4/Refmac monomer CIF library.
///
/// # Errors
///
/// Returns an error for duplicate components, unresolved atoms, or malformed
/// bond, angle, torsion, plane, or chirality restraints.
pub fn lower_monomer_library(document: &Document) -> Result<MonomerLibrary, RestraintError> {
    let mut components = BTreeMap::new();
    for block in document.blocks() {
        if block.category("chem_comp_atom").is_none() {
            continue;
        }
        let component = lower_component(block)?;
        if components.insert(component.id.clone(), component).is_some() {
            return Err(RestraintError::ComponentId);
        }
    }
    if components.is_empty() {
        return Err(RestraintError::MissingComponents);
    }
    Ok(MonomerLibrary { components })
}

/// Parses a CIF monomer library through the shared CIF parser.
///
/// # Errors
///
/// Returns parser diagnostics separately from restraint semantic errors.
pub fn read_monomer_library(
    input: &InputBuffer,
) -> Result<(MonomerLibrary, Vec<pdbiox_core::Diagnostic>), MonomerLibraryReadError> {
    let (document, findings) =
        pdbiox_cif::parse(input).map_err(MonomerLibraryReadError::Findings)?;
    lower_monomer_library(&document)
        .map(|library| (library, findings))
        .map_err(MonomerLibraryReadError::Restraints)
}

/// Monomer-library syntax or semantic failure.
#[derive(Debug)]
#[non_exhaustive]
pub enum MonomerLibraryReadError {
    /// CIF parser diagnostics.
    Findings(Vec<pdbiox_core::Diagnostic>),
    /// Restraint semantic error.
    Restraints(RestraintError),
}
impl std::fmt::Display for MonomerLibraryReadError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Findings(values) => write!(
                formatter,
                "monomer CIF parsing failed with {} findings",
                values.len()
            ),
            Self::Restraints(error) => error.fmt(formatter),
        }
    }
}
impl std::error::Error for MonomerLibraryReadError {}

fn lower_component(block: &DataBlock) -> Result<MonomerRestraints, RestraintError> {
    let id = match block
        .category("chem_comp")
        .and_then(|category| category.identifier("id", 0))
    {
        Some(value) => value.into_owned(),
        None => block.name().trim_start_matches("comp_").to_owned(),
    };
    if id.is_empty() {
        return Err(RestraintError::ComponentId);
    }
    let atom_category = block
        .category("chem_comp_atom")
        .ok_or(RestraintError::MissingComponents)?;
    let atoms: Vec<Box<str>> = (0..atom_category.row_count())
        .map(|row| identifier(atom_category, "atom_id", row))
        .collect::<Result<_, _>>()?;
    let known: BTreeSet<&str> = atoms.iter().map(Box::as_ref).collect();
    Ok(MonomerRestraints {
        id: id.into(),
        bonds: lower_bonds(block, &known)?,
        angles: lower_angles(block, &known)?,
        torsions: lower_torsions(block, &known)?,
        planes: lower_planes(block, &known)?,
        chirals: lower_chirals(block, &known)?,
        atoms,
    })
}

fn lower_bonds(
    block: &DataBlock,
    known: &BTreeSet<&str>,
) -> Result<Vec<BondRestraint>, RestraintError> {
    let Some(category) = block.category("chem_comp_bond") else {
        return Ok(Vec::new());
    };
    (0..category.row_count())
        .map(|row| {
            let atoms = [
                identifier(category, "atom_id_1", row)?,
                identifier(category, "atom_id_2", row)?,
            ];
            check_atoms(&atoms, known)?;
            Ok(BondRestraint {
                atoms,
                target: positive_number(category, "value_dist", row)?,
                sigma: positive_number(category, "value_dist_esd", row)?,
                kind: category
                    .identifier("type", row)
                    .map(|value| value.into_owned().into()),
            })
        })
        .collect()
}
fn lower_angles(
    block: &DataBlock,
    known: &BTreeSet<&str>,
) -> Result<Vec<AngleRestraint>, RestraintError> {
    let Some(category) = block.category("chem_comp_angle") else {
        return Ok(Vec::new());
    };
    (0..category.row_count())
        .map(|row| {
            let atoms = [
                identifier(category, "atom_id_1", row)?,
                identifier(category, "atom_id_2", row)?,
                identifier(category, "atom_id_3", row)?,
            ];
            check_atoms(&atoms, known)?;
            Ok(AngleRestraint {
                atoms,
                target: finite_number(category, "value_angle", row)?,
                sigma: positive_number(category, "value_angle_esd", row)?,
            })
        })
        .collect()
}
fn lower_torsions(
    block: &DataBlock,
    known: &BTreeSet<&str>,
) -> Result<Vec<TorsionRestraint>, RestraintError> {
    let Some(category) = block.category("chem_comp_tor") else {
        return Ok(Vec::new());
    };
    (0..category.row_count())
        .map(|row| {
            let atoms = [
                identifier(category, "atom_id_1", row)?,
                identifier(category, "atom_id_2", row)?,
                identifier(category, "atom_id_3", row)?,
                identifier(category, "atom_id_4", row)?,
            ];
            check_atoms(&atoms, known)?;
            let period = category
                .value("period", row)
                .and_then(CifValue::as_integer)
                .and_then(|value| u32::try_from(value).ok())
                .filter(|value| *value > 0)
                .ok_or(RestraintError::InvalidValue)?;
            Ok(TorsionRestraint {
                id: category
                    .identifier("id", row)
                    .or_else(|| category.identifier("label", row))
                    .map(|value| value.into_owned().into()),
                atoms,
                target: finite_number(category, "value_angle", row)?,
                sigma: positive_number(category, "value_angle_esd", row)?,
                period,
            })
        })
        .collect()
}
fn lower_planes(
    block: &DataBlock,
    known: &BTreeSet<&str>,
) -> Result<Vec<PlaneRestraint>, RestraintError> {
    let Some(category) = block
        .category("chem_comp_plane_atom")
        .or_else(|| block.category("chem_comp_plane"))
    else {
        return Ok(Vec::new());
    };
    let id_item = if category.column("plane_id").is_some() {
        "plane_id"
    } else {
        "id"
    };
    let sigma_item = if category.column("dist_esd").is_some() {
        "dist_esd"
    } else {
        "esd"
    };
    let mut groups: BTreeMap<Box<str>, Vec<PlaneAtomRestraint>> = BTreeMap::new();
    for row in 0..category.row_count() {
        let id = identifier(category, id_item, row)?;
        let atom = identifier(category, "atom_id", row)?;
        check_atoms(std::slice::from_ref(&atom), known)?;
        groups.entry(id).or_default().push(PlaneAtomRestraint {
            atom,
            sigma: positive_number(category, sigma_item, row)?,
        });
    }
    Ok(groups
        .into_iter()
        .map(|(id, atoms)| PlaneRestraint { id, atoms })
        .collect())
}
fn lower_chirals(
    block: &DataBlock,
    known: &BTreeSet<&str>,
) -> Result<Vec<ChiralRestraint>, RestraintError> {
    let Some(category) = block.category("chem_comp_chir") else {
        return Ok(Vec::new());
    };
    (0..category.row_count())
        .map(|row| {
            let modern = category.column("atom_id_centre").is_some();
            let atoms = if modern {
                [
                    identifier(category, "atom_id_centre", row)?,
                    identifier(category, "atom_id_1", row)?,
                    identifier(category, "atom_id_2", row)?,
                    identifier(category, "atom_id_3", row)?,
                ]
            } else {
                [
                    identifier(category, "atom_id_1", row)?,
                    identifier(category, "atom_id_2", row)?,
                    identifier(category, "atom_id_3", row)?,
                    identifier(category, "atom_id_4", row)?,
                ]
            };
            check_atoms(&atoms, known)?;
            let sign = match category
                .text("volume_sign", row)
                .map(str::to_ascii_lowercase)
                .as_deref()
            {
                Some("positive" | "positiv") => ChiralVolumeSign::Positive,
                Some("negative" | "negativ") => ChiralVolumeSign::Negative,
                Some("both") => ChiralVolumeSign::Both,
                _ => return Err(RestraintError::InvalidValue),
            };
            Ok(ChiralRestraint {
                id: category
                    .identifier("id", row)
                    .map(|value| value.into_owned().into()),
                atoms,
                sign,
            })
        })
        .collect()
}

fn identifier(
    category: &pdbiox_cif::Category,
    item: &str,
    row: usize,
) -> Result<Box<str>, RestraintError> {
    category
        .identifier(item, row)
        .map(|value| value.into_owned().into())
        .ok_or(RestraintError::InvalidValue)
}
fn finite_number(
    category: &pdbiox_cif::Category,
    item: &str,
    row: usize,
) -> Result<f64, RestraintError> {
    category
        .value(item, row)
        .and_then(CifValue::as_float)
        .filter(|value| value.is_finite())
        .ok_or(RestraintError::InvalidValue)
}
fn positive_number(
    category: &pdbiox_cif::Category,
    item: &str,
    row: usize,
) -> Result<f64, RestraintError> {
    finite_number(category, item, row).and_then(|value| {
        if value > 0.0 {
            Ok(value)
        } else {
            Err(RestraintError::InvalidValue)
        }
    })
}
fn check_atoms(atoms: &[Box<str>], known: &BTreeSet<&str>) -> Result<(), RestraintError> {
    if atoms.iter().all(|atom| known.contains(atom.as_ref())) {
        Ok(())
    } else {
        Err(RestraintError::AtomReference)
    }
}

#[cfg(test)]
#[path = "restraints_tests.rs"]
mod tests;
