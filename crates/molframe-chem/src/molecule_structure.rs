//! Conversion between small-molecule graphs and structures.
//!
//! A molecule becomes one model holding one chain `A` that holds one residue,
//! the way a ligand file's `HETATM` block reads. Every bond keeps the order the
//! file declared and the provenance `File`.

use crate::{Mol2Record, MolAtom, MolBond, MolRecord, Molecule};
use molframe_cif::SmallCifStructure;
use molframe_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::chunk::{AtomRecord, ChunkBuilder};
use molframe_core::column::Presence;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::element::Element;
use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::optional::{OptionalI32, OptionalSymbol};
use molframe_core::secondary::{SecondarySource, SecondaryStructure};
use molframe_core::structure::{CoordinateStore, Structure, StructureData, UnitCell};
use molframe_core::symbol::{AltId, SymbolId};
use molframe_core::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};

/// The residue name used when a molecule carries no usable name.
const UNNAMED_RESIDUE: &str = "UNL";

/// One atom to place, with the label and charge its source declared.
struct Site {
    atom: MolAtom,
    name: String,
    formal_charge: Option<i8>,
}

fn consistency(message: &str) -> Diagnostic {
    Diagnostic::new(Code::E3001).with_message(message.to_owned())
}

fn full<T>(_: T) -> Diagnostic {
    consistency("the identifier dictionary is full")
}

fn residue_name(name: &str) -> String {
    let text: String = name
        .trim()
        .chars()
        .filter(char::is_ascii_alphanumeric)
        .take(3)
        .collect::<String>()
        .to_ascii_uppercase();
    if text.is_empty() {
        UNNAMED_RESIDUE.to_owned()
    } else {
        text
    }
}

/// The unique per-element label `C1`, `C2`, ... for atoms with no name.
fn generated_names(molecule: &Molecule) -> Vec<String> {
    let mut counts = std::collections::HashMap::new();
    molecule
        .atoms
        .iter()
        .map(|atom| {
            let symbol = atom.element.symbol();
            let count = counts.entry(symbol).or_insert(0_usize);
            *count += 1;
            format!("{symbol}{count}")
        })
        .collect()
}

fn bond_order(order: u8) -> BondOrder {
    match order {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        _ => BondOrder::Unknown,
    }
}

fn add_atoms(data: &mut StructureData, sites: &[Site]) -> Result<(), Diagnostic> {
    let mut builder = ChunkBuilder::new();
    for (position, site) in sites.iter().enumerate() {
        let atom_name = data.dictionary.intern(&site.name).map_err(full)?;
        let serial = u32::try_from(position + 1)
            .map_err(|_| consistency("the molecule exceeds the supported index range"))?;
        builder.push(AtomRecord {
            position: Some(site.atom.position),
            element: site.atom.element,
            atom_name,
            auth_atom_name: OptionalSymbol::some(atom_name),
            alternate_component_id: OptionalSymbol::NONE,
            alt_id: AltId::BLANK,
            residue: ResidueIndex::new(0),
            occupancy: (1.0, Presence::Unknown),
            b_factor: (0.0, Presence::Unknown),
            formal_charge: match site.formal_charge {
                Some(charge) => (charge, Presence::Present),
                None => (0, Presence::Unknown),
            },
            atom_site_id: serial,
        });
    }
    let (chunks, coordinates) = builder.finish();
    data.chunks = chunks.into();
    data.coords = CoordinateStore::Single(coordinates);
    Ok(())
}

fn add_hierarchy(
    data: &mut StructureData,
    component: SymbolId,
    atoms: u32,
) -> Result<(), Diagnostic> {
    let chain = data.dictionary.intern("A").map_err(full)?;
    let entity_id = data.dictionary.intern("1").map_err(full)?;
    let failed = |error: &dyn std::fmt::Display| consistency(&error.to_string());
    data.topology
        .residues
        .push(
            ResidueRecord {
                label_comp_id: component,
                auth_comp_id: OptionalSymbol::some(component),
                label_seq_id: OptionalI32::some(1),
                auth_seq_id: OptionalI32::some(1),
                ins_code: OptionalSymbol::NONE,
                het: true,
            },
            0..atoms,
        )
        .map_err(|error| failed(&error))?;
    let entity = data
        .topology
        .entities
        .push(entity_id, EntityKind::Unknown, OptionalSymbol::NONE, &[])
        .map_err(|error| failed(&error))?;
    data.topology
        .chains
        .push(
            ChainRecord {
                label_asym_id: chain,
                auth_asym_id: OptionalSymbol::some(chain),
                entity,
                polymer_kind: PolymerKind::None,
            },
            0..1,
        )
        .map_err(|error| failed(&error))?;
    data.topology
        .models
        .push(1, 0..1)
        .map_err(|error| failed(&error))?;
    Ok(())
}

fn build(
    name: &str,
    sites: &[Site],
    bonds: &[(usize, usize, BondOrder)],
) -> Result<Structure, Diagnostic> {
    let mut data = StructureData::empty();
    let atoms = u32::try_from(sites.len())
        .map_err(|_| consistency("the molecule exceeds the supported index range"))?;
    let component = data.dictionary.intern(&residue_name(name)).map_err(full)?;
    add_atoms(&mut data, sites)?;
    add_hierarchy(&mut data, component, atoms)?;

    let mut table = BondTableBuilder::new();
    for &(first, second, order) in bonds {
        let (Ok(first), Ok(second)) = (u32::try_from(first), u32::try_from(second)) else {
            return Err(consistency("a bond endpoint exceeds the supported range"));
        };
        if first >= atoms || second >= atoms || first == second {
            return Err(consistency("a bond names an atom the molecule lacks"));
        }
        table.push(BondRecord {
            atom_a: AtomIndex::new(first),
            atom_b: AtomIndex::new(second),
            order,
            provenance: BondProvenance::File,
        });
    }
    data.bonds = table.finish();
    data.secondary_structure = vec![SecondaryStructure::Unknown].into();
    data.secondary_source = vec![SecondarySource::None].into();
    Ok(Structure::new(data))
}

fn plain_sites(molecule: &Molecule) -> Vec<Site> {
    generated_names(molecule)
        .into_iter()
        .zip(&molecule.atoms)
        .map(|(name, atom)| Site {
            atom: *atom,
            name,
            formal_charge: None,
        })
        .collect()
}

fn plain_bonds(bonds: &[MolBond]) -> Vec<(usize, usize, BondOrder)> {
    bonds
        .iter()
        .map(|bond| (bond.first, bond.second, bond_order(bond.order)))
        .collect()
}

/// Lowers a molecule to a one-residue structure named `name`.
///
/// Atoms are labelled by element and running count (`C1`, `C2`, `N1`, ...).
///
/// # Errors
///
/// Returns a diagnostic when a bond names an absent atom or the molecule
/// exceeds the structure index space.
pub fn molecule_to_structure(name: &str, molecule: &Molecule) -> Result<Structure, Diagnostic> {
    build(name, &plain_sites(molecule), &plain_bonds(&molecule.bonds))
}

/// Lowers a MOL block or SDF record, keeping its declared formal charges.
///
/// # Errors
///
/// As [`molecule_to_structure`].
pub fn mol_record_to_structure(record: &MolRecord) -> Result<Structure, Diagnostic> {
    let mut sites = plain_sites(&record.molecule);
    for (site, metadata) in sites.iter_mut().zip(&record.atom_metadata) {
        site.formal_charge = metadata.formal_charge;
    }
    build(&record.name, &sites, &plain_bonds(&record.molecule.bonds))
}

/// Lowers a MOL2 record, keeping its atom names and bond types.
///
/// # Errors
///
/// As [`molecule_to_structure`], and when the record's metadata does not match
/// its graph.
pub fn mol2_record_to_structure(record: &Mol2Record) -> Result<Structure, Diagnostic> {
    if record.atom_metadata.len() != record.molecule.atoms.len()
        || record.bond_metadata.len() != record.molecule.bonds.len()
    {
        return Err(consistency("MOL2 metadata does not match its graph"));
    }
    let sites: Vec<Site> = record
        .molecule
        .atoms
        .iter()
        .zip(&record.atom_metadata)
        .map(|(atom, metadata)| Site {
            atom: *atom,
            name: metadata.name.to_string(),
            formal_charge: None,
        })
        .collect();
    let bonds: Vec<(usize, usize, BondOrder)> = record
        .molecule
        .bonds
        .iter()
        .zip(&record.bond_metadata)
        .map(|(bond, metadata)| {
            let order = match &*metadata.bond_type {
                "1" | "am" => BondOrder::Single,
                "2" => BondOrder::Double,
                "3" => BondOrder::Triple,
                "ar" => BondOrder::Aromatic,
                _ => BondOrder::Unknown,
            };
            (bond.first, bond.second, order)
        })
        .collect();
    build(&record.name, &sites, &bonds)
}

/// Extracts a structure's atoms and bonds as one molecule graph.
///
/// # Errors
///
/// Returns a diagnostic when an atom has no coordinates or a bond's order is
/// unknown, since neither can be written without inventing a value.
pub fn structure_to_molecule(structure: &Structure) -> Result<Molecule, Diagnostic> {
    let mut atoms = Vec::new();
    for atom in structure.data().atoms() {
        let position = atom
            .position()
            .ok_or_else(|| consistency("an atom has no recorded coordinates"))?;
        let element = atom
            .element()
            .ok_or_else(|| consistency("an atom has no element"))?;
        atoms.push(MolAtom { element, position });
    }
    let mut bonds = Vec::with_capacity(structure.data().bonds.len());
    for record in structure.data().bonds.iter() {
        let order = match record.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            BondOrder::Quadruple | BondOrder::Polymeric | BondOrder::Unknown => {
                return Err(consistency("a bond has no representable order"));
            }
        };
        let (Ok(first), Ok(second)) = (
            usize::try_from(record.atom_a.get()),
            usize::try_from(record.atom_b.get()),
        ) else {
            return Err(consistency("a bond endpoint exceeds the supported range"));
        };
        bonds.push(MolBond {
            first,
            second,
            order,
        });
    }
    Ok(Molecule { atoms, bonds })
}

/// The Cartesian position of a fractional site in a cell (International
/// Tables orientation: `a` along x, `b` in the xy plane).
fn orthogonalise(cell: &UnitCell, fractional: [f64; 3]) -> Option<[f64; 3]> {
    let [alpha, beta, gamma] = cell.angles.map(f64::to_radians);
    let (cos_alpha, cos_beta, cos_gamma) = (alpha.cos(), beta.cos(), gamma.cos());
    let sin_gamma = gamma.sin();
    let volume_term = 1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma
        + 2.0 * cos_alpha * cos_beta * cos_gamma;
    if sin_gamma.abs() < 1e-12 || volume_term <= 0.0 {
        return None;
    }
    // Columns of the orthogonalisation matrix: one cell edge each.
    let [length_a, length_b, length_c] = cell.lengths;
    let edge_a = [length_a, 0.0, 0.0];
    let edge_b = [length_b * cos_gamma, length_b * sin_gamma, 0.0];
    let edge_c = [
        length_c * cos_beta,
        length_c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma,
        length_c * volume_term.sqrt() / sin_gamma,
    ];
    let mut position = [0.0; 3];
    for (axis, value) in position.iter_mut().enumerate() {
        *value = edge_a[axis] * fractional[0]
            + edge_b[axis] * fractional[1]
            + edge_c[axis] * fractional[2];
    }
    Some(position)
}

/// The element a deposited type symbol names, ignoring a trailing oxidation
/// state such as `Fe3+` or `O2-`.
fn element_of(type_symbol: &str) -> Option<Element> {
    let letters: String = type_symbol
        .chars()
        .take_while(char::is_ascii_alphabetic)
        .collect();
    let mut candidates = Vec::new();
    if letters.len() >= 2 {
        candidates.push(letters[..2].to_owned());
    }
    candidates.extend(letters.get(..1).map(str::to_owned));
    candidates
        .iter()
        .find_map(|symbol| Element::from_symbol(symbol))
}

/// Lowers a core-CIF crystal model to a one-residue structure.
///
/// Cartesian coordinates are used where the file records them; otherwise
/// fractional ones are orthogonalised with the unit cell. Bonds keep an
/// unknown order, because core CIF records distances rather than orders.
///
/// # Errors
///
/// Returns a diagnostic when an atom's element or position cannot be
/// determined, or the cell is degenerate.
pub fn small_cif_to_structure(model: &SmallCifStructure) -> Result<Structure, Diagnostic> {
    let mut sites = Vec::with_capacity(model.atoms.len());
    for atom in &model.atoms {
        let element = element_of(&atom.type_symbol)
            .ok_or_else(|| consistency("a site names no known element"))?;
        let cartesian = match (atom.cartesian, atom.fractional, model.cell.as_ref()) {
            (Some(cartesian), _, _) => cartesian,
            (None, Some(fractional), Some(cell)) => orthogonalise(cell, fractional)
                .ok_or_else(|| consistency("the unit cell is degenerate"))?,
            _ => return Err(consistency("a site has no usable coordinates")),
        };
        #[allow(clippy::cast_possible_truncation)]
        let position = cartesian.map(|value| value as f32);
        sites.push(Site {
            atom: MolAtom { element, position },
            name: atom.label.to_string(),
            formal_charge: None,
        });
    }
    let bonds: Vec<(usize, usize, BondOrder)> = model
        .bonds
        .iter()
        .map(|bond| (bond.first, bond.second, BondOrder::Unknown))
        .collect();
    build(&model.name, &sites, &bonds)
}

#[cfg(test)]
#[path = "molecule_structure_tests.rs"]
mod tests;
