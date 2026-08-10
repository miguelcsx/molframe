//! Molecular-topology tables and molecule-block expansion.

use super::constants::{
    DIMENSIONS, F_SETTLE, FUNCTION_COUNT, function_absent, interaction_arity, is_bond_function,
    remap_function,
};
use super::error::TprError;
use super::model::{TprAtom, TprBond, TprHeader, TprResidue, TprTopology};
use super::parameters::skip_parameter_table;
use super::xdr::Decoder;

#[derive(Clone, Debug)]
struct AtomTemplate {
    name: String,
    atom_type: String,
    residue: usize,
    mass: f64,
    charge: f64,
    atomic_number: Option<u8>,
}

#[derive(Clone, Debug)]
struct MoleculeTemplate {
    name: String,
    residue_names: Vec<String>,
    atoms: Vec<AtomTemplate>,
    bonds: Vec<TprBond>,
}

#[derive(Clone, Copy, Debug)]
struct MoleculeBlock {
    template: usize,
    count: usize,
    atom_count: usize,
}

pub(super) fn parse_topology(
    decoder: &mut Decoder<'_>,
    header: &TprHeader,
) -> Result<TprTopology, TprError> {
    let symbols = read_symbols(decoder)?;
    read_symbol(decoder, &symbols, "system name")?;
    let _atom_type_count = read_force_field(decoder, header.format_version)?;
    let template_count = decoder.count("molecule type count")?;
    let mut templates = Vec::with_capacity(template_count);
    for _ in 0..template_count {
        templates.push(read_molecule_template(
            decoder,
            &symbols,
            header.format_version,
        )?);
    }
    let block_count = decoder.count("molecule block count")?;
    let mut blocks = Vec::with_capacity(block_count);
    for _ in 0..block_count {
        blocks.push(read_molecule_block(decoder)?);
    }
    expand_topology(header, &templates, &blocks)
}

fn read_symbols(decoder: &mut Decoder<'_>) -> Result<Vec<String>, TprError> {
    let count = decoder.count("symbol count")?;
    (0..count).map(|_| decoder.string("symbol")).collect()
}

fn read_symbol<'a>(
    decoder: &mut Decoder<'_>,
    symbols: &'a [String],
    field: &'static str,
) -> Result<&'a str, TprError> {
    let index = decoder.count(field)?;
    symbols
        .get(index)
        .map(String::as_str)
        .ok_or(TprError::InvalidIndex {
            field,
            index,
            length: symbols.len(),
        })
}

fn read_force_field(decoder: &mut Decoder<'_>, version: i32) -> Result<usize, TprError> {
    let atom_type_count = decoder.count("force-field atom type count")?;
    let parameter_count = decoder.count("force-field parameter count")?;
    let mut functions = Vec::with_capacity(parameter_count);
    for _ in 0..parameter_count {
        functions.push(remap_function(version, decoder.i32()?));
    }
    if version >= 66 {
        decoder.double()?;
    }
    decoder.real()?;
    skip_parameter_table(decoder, &functions, version)?;
    Ok(atom_type_count)
}

fn read_molecule_template(
    decoder: &mut Decoder<'_>,
    symbols: &[String],
    version: i32,
) -> Result<MoleculeTemplate, TprError> {
    let name = read_symbol(decoder, symbols, "molecule type name")?.to_owned();
    let atom_count = decoder.count("molecule atom count")?;
    let residue_count = decoder.count("molecule residue count")?;
    let mut raw_atoms = Vec::with_capacity(atom_count);
    for _ in 0..atom_count {
        raw_atoms.push(read_atom_template(decoder)?);
    }
    let atom_names = read_symbol_array(decoder, symbols, atom_count, "atom name")?;
    let atom_types = read_symbol_array(decoder, symbols, atom_count, "atom type")?;
    let _state_b_types = read_symbol_array(decoder, symbols, atom_count, "state B atom type")?;
    let residue_names = read_residue_names(decoder, symbols, residue_count, version)?;
    let atoms = assemble_atoms(raw_atoms, atom_names, atom_types, residue_names.len())?;
    let bonds = read_interactions(decoder, version)?;
    skip_block(decoder)?;
    skip_blocka(decoder)?;
    Ok(MoleculeTemplate {
        name,
        residue_names,
        atoms,
        bonds,
    })
}

#[derive(Clone, Copy, Debug)]
struct RawAtom {
    residue: usize,
    mass: f64,
    charge: f64,
    atomic_number: Option<u8>,
}

fn read_atom_template(decoder: &mut Decoder<'_>) -> Result<RawAtom, TprError> {
    let mass = decoder.real()?;
    let charge = decoder.real()?;
    decoder.skip_reals(2)?;
    decoder.ushort()?;
    decoder.ushort()?;
    decoder.i32()?;
    let residue = decoder.count("atom residue index")?;
    let atomic_number = match decoder.i32()? {
        value @ 1..=118 => Some(u8::try_from(value).map_err(|_| TprError::SizeOverflow)?),
        _ => None,
    };
    Ok(RawAtom {
        residue,
        mass,
        charge,
        atomic_number,
    })
}

fn read_symbol_array(
    decoder: &mut Decoder<'_>,
    symbols: &[String],
    count: usize,
    field: &'static str,
) -> Result<Vec<String>, TprError> {
    (0..count)
        .map(|_| read_symbol(decoder, symbols, field).map(str::to_owned))
        .collect()
}

fn read_residue_names(
    decoder: &mut Decoder<'_>,
    symbols: &[String],
    count: usize,
    version: i32,
) -> Result<Vec<String>, TprError> {
    let mut names = Vec::with_capacity(count);
    for _ in 0..count {
        names.push(read_symbol(decoder, symbols, "residue name")?.to_owned());
        if version >= 63 {
            decoder.i32()?;
            decoder.uchar()?;
        }
    }
    Ok(names)
}

fn assemble_atoms(
    raw_atoms: Vec<RawAtom>,
    names: Vec<String>,
    atom_types: Vec<String>,
    residue_count: usize,
) -> Result<Vec<AtomTemplate>, TprError> {
    raw_atoms
        .into_iter()
        .zip(names)
        .zip(atom_types)
        .map(|((raw, name), atom_type)| {
            if raw.residue >= residue_count {
                return Err(TprError::InvalidIndex {
                    field: "atom residue",
                    index: raw.residue,
                    length: residue_count,
                });
            }
            Ok(AtomTemplate {
                name,
                atom_type,
                residue: raw.residue,
                mass: raw.mass,
                charge: raw.charge,
                atomic_number: raw.atomic_number,
            })
        })
        .collect()
}

fn read_interactions(decoder: &mut Decoder<'_>, version: i32) -> Result<Vec<TprBond>, TprError> {
    let mut bonds = Vec::new();
    for function in 0..FUNCTION_COUNT {
        if function_absent(version, function) {
            continue;
        }
        let value_count = decoder.count("interaction value count")?;
        let mut values = Vec::with_capacity(value_count);
        for _ in 0..value_count {
            values.push(decoder.i32()?);
        }
        append_bonds(&mut bonds, function, &values)?;
    }
    bonds.sort_unstable_by_key(|bond| (bond.atom_a, bond.atom_b));
    bonds.dedup();
    Ok(bonds)
}

fn append_bonds(bonds: &mut Vec<TprBond>, function: i32, values: &[i32]) -> Result<(), TprError> {
    if function == F_SETTLE {
        return append_settle_bonds(bonds, values);
    }
    if !is_bond_function(function) {
        return Ok(());
    }
    let Some(arity) = interaction_arity(function) else {
        return Err(TprError::UnsupportedFunction(function));
    };
    for record in interaction_records(values, arity, "bond interaction")? {
        bonds.push(ordered_bond(atom_index(record[0])?, atom_index(record[1])?));
    }
    Ok(())
}

fn append_settle_bonds(bonds: &mut Vec<TprBond>, values: &[i32]) -> Result<(), TprError> {
    if values.len() == 2 {
        let oxygen = atom_index(values[1])?;
        bonds.push(ordered_bond(oxygen, oxygen + 1));
        bonds.push(ordered_bond(oxygen, oxygen + 2));
        return Ok(());
    }
    for record in interaction_records(values, 3, "SETTLE interaction")? {
        let oxygen = atom_index(record[0])?;
        bonds.push(ordered_bond(oxygen, atom_index(record[1])?));
        bonds.push(ordered_bond(oxygen, atom_index(record[2])?));
    }
    Ok(())
}

fn interaction_records<'a>(
    values: &'a [i32],
    arity: usize,
    field: &'static str,
) -> Result<impl Iterator<Item = &'a [i32]>, TprError> {
    let width = arity + 1;
    if !values.len().is_multiple_of(width) {
        let value = i64::try_from(values.len()).map_err(|_| TprError::SizeOverflow)?;
        return Err(TprError::InvalidCount {
            field,
            value,
            offset: 0,
        });
    }
    Ok(values.chunks_exact(width).map(|record| &record[1..]))
}

fn atom_index(value: i32) -> Result<usize, TprError> {
    usize::try_from(value).map_err(|_| TprError::InvalidCount {
        field: "interaction atom index",
        value: i64::from(value),
        offset: 0,
    })
}

const fn ordered_bond(atom_a: usize, atom_b: usize) -> TprBond {
    if atom_a <= atom_b {
        TprBond { atom_a, atom_b }
    } else {
        TprBond {
            atom_a: atom_b,
            atom_b: atom_a,
        }
    }
}

fn skip_block(decoder: &mut Decoder<'_>) -> Result<(), TprError> {
    let count = decoder.count("charge-group block count")?;
    skip_ints(decoder, count.checked_add(1).ok_or(TprError::SizeOverflow)?)
}

fn skip_blocka(decoder: &mut Decoder<'_>) -> Result<(), TprError> {
    let block_count = decoder.count("exclusion block count")?;
    let value_count = decoder.count("exclusion value count")?;
    skip_ints(
        decoder,
        block_count.checked_add(1).ok_or(TprError::SizeOverflow)?,
    )?;
    skip_ints(decoder, value_count)
}

fn skip_ints(decoder: &mut Decoder<'_>, count: usize) -> Result<(), TprError> {
    for _ in 0..count {
        decoder.i32()?;
    }
    Ok(())
}

fn read_molecule_block(decoder: &mut Decoder<'_>) -> Result<MoleculeBlock, TprError> {
    let template = decoder.count("molecule block type")?;
    let count = decoder.count("molecule block count")?;
    let atom_count = decoder.count("molecule block atom count")?;
    let position_count = decoder.count("position-restraint A count")?;
    decoder.skip_reals(
        position_count
            .checked_mul(DIMENSIONS)
            .ok_or(TprError::SizeOverflow)?,
    )?;
    let alternate_position_count = decoder.count("position-restraint B count")?;
    decoder.skip_reals(
        alternate_position_count
            .checked_mul(DIMENSIONS)
            .ok_or(TprError::SizeOverflow)?,
    )?;
    Ok(MoleculeBlock {
        template,
        count,
        atom_count,
    })
}

fn expand_topology(
    header: &TprHeader,
    templates: &[MoleculeTemplate],
    blocks: &[MoleculeBlock],
) -> Result<TprTopology, TprError> {
    let mut topology = TprTopology {
        header: header.clone(),
        residues: Vec::new(),
        atoms: Vec::with_capacity(header.atom_count),
        bonds: Vec::new(),
    };
    let mut molecule = 0;
    for block in blocks {
        let template = templates
            .get(block.template)
            .ok_or(TprError::InvalidIndex {
                field: "molecule template",
                index: block.template,
                length: templates.len(),
            })?;
        if template.atoms.len() != block.atom_count {
            let value = i64::try_from(block.atom_count).map_err(|_| TprError::SizeOverflow)?;
            return Err(TprError::InvalidCount {
                field: "molecule block atom count",
                value,
                offset: 0,
            });
        }
        for _ in 0..block.count {
            append_molecule(&mut topology, template, molecule)?;
            molecule += 1;
        }
    }
    if topology.atoms.len() != header.atom_count {
        let value = i64::try_from(topology.atoms.len()).map_err(|_| TprError::SizeOverflow)?;
        return Err(TprError::InvalidCount {
            field: "expanded atom count",
            value,
            offset: 0,
        });
    }
    topology
        .bonds
        .sort_unstable_by_key(|bond| (bond.atom_a, bond.atom_b));
    topology.bonds.dedup();
    Ok(topology)
}

fn append_molecule(
    topology: &mut TprTopology,
    template: &MoleculeTemplate,
    molecule: usize,
) -> Result<(), TprError> {
    let residue_offset = topology.residues.len();
    let atom_offset = topology.atoms.len();
    for (position, name) in template.residue_names.iter().enumerate() {
        let number =
            i32::try_from(residue_offset + position + 1).map_err(|_| TprError::SizeOverflow)?;
        topology.residues.push(TprResidue {
            name: name.clone(),
            number,
            molecule,
            molecule_type: template.name.clone(),
        });
    }
    for atom in &template.atoms {
        topology.atoms.push(TprAtom {
            index: topology.atoms.len(),
            name: atom.name.clone(),
            atom_type: atom.atom_type.clone(),
            residue: residue_offset + atom.residue,
            mass: atom.mass,
            charge: atom.charge,
            atomic_number: atom.atomic_number,
        });
    }
    for bond in &template.bonds {
        topology.bonds.push(ordered_bond(
            atom_offset + bond.atom_a,
            atom_offset + bond.atom_b,
        ));
    }
    Ok(())
}
