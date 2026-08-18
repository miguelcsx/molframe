/// Returns whether two structures have identical normalized topology.
///
/// The comparison covers declared hierarchy sizes, bonds, chain metadata,
/// residue identity fields, and atom identity fields. Coordinate-dependent
/// properties are intentionally excluded.
fn same_topology(left: &Structure, right: &Structure) -> bool {
    if left.atom_count() != right.atom_count()
        || left.residue_count() != right.residue_count()
        || left.chain_count() != right.chain_count()
    {
        return false;
    }

    let left_data = left.data();
    let right_data = right.data();

    left_data.bonds.iter().eq(right_data.bonds.iter())
        && same_chain_topology(left_data, right_data)
        && same_residue_topology(left_data, right_data)
        && same_atom_topology(left_data, right_data)
}

/// Compares the topology-relevant metadata of every chain in two structures.
///
/// Chain labels, author labels, and polymer classifications must match in
/// iteration order.
fn same_chain_topology(left: &StructureData, right: &StructureData) -> bool {
    left.chains().zip(right.chains()).all(|(left, right)| {
        left.label() == right.label()
            && left.auth_label() == right.auth_label()
            && left.polymer_kind() == right.polymer_kind()
    })
}

/// Compares the topology-relevant metadata of every residue in two structures.
///
/// Residue names, sequence identifiers, insertion codes, and hetero flags must
/// match in iteration order.
fn same_residue_topology(left: &StructureData, right: &StructureData) -> bool {
    left.residues().zip(right.residues()).all(|(left, right)| {
        left.name() == right.name()
            && left.auth_seq_id() == right.auth_seq_id()
            && left.label_seq_id() == right.label_seq_id()
            && left.ins_code() == right.ins_code()
            && left.is_het() == right.is_het()
    })
}

/// Compares the topology-relevant metadata of every atom in two structures.
///
/// Atom names, elements, alternate labels, and formal charges must match in
/// iteration order.
fn same_atom_topology(left: &StructureData, right: &StructureData) -> bool {
    left.atoms().zip(right.atoms()).all(|(left, right)| {
        left.name() == right.name()
            && left.element() == right.element()
            && left.alt_label() == right.alt_label()
            && left.formal_charge() == right.formal_charge()
    })
}

/// Builds one normalized atom record from an MMTF group entry.
///
/// Resolves interned identifiers, coordinates, alternate locations, optional
/// experimental values, formal charge, and the source atom-site identifier.
fn atom_record(
    data: &mut StructureData,
    decoded: &Decoded,
    group: &Group,
    local: usize,
    atom: usize,
    residue: ResidueIndex,
) -> Result<AtomRecord, Diagnostic> {
    let atom_name = intern(
        data,
        group
            .atom_name_list
            .get(local)
            .ok_or_else(|| schema_error("group atom name missing"))?,
    )?;

    let alt_symbol = char_symbol(
        data,
        decoded
            .alt_loc
            .as_ref()
            .and_then(|values| values.get(atom))
            .copied(),
    )?;

    let alt_id = match alt_symbol {
        Some(symbol) => AltId::labelled(symbol)
            .ok_or_else(|| schema_error("alternate-location identifier exceeds its encoding"))?,
        None => AltId::BLANK,
    };

    let charge = group
        .formal_charge_list
        .get(local)
        .copied()
        .ok_or_else(|| schema_error("group formal charge missing"))?;

    let formal_charge =
        i8::try_from(charge).map_err(|_| schema_error("formal charge is outside i8"))?;

    Ok(AtomRecord {
        position: Some(position(decoded, atom)?),
        element: group_element(group, local),
        atom_name,
        auth_atom_name: OptionalSymbol::NONE,
        alternate_component_id: OptionalSymbol::NONE,
        alt_id,
        residue,
        occupancy: optional_float(decoded.occupancy.as_ref(), atom, 1.0),
        b_factor: optional_float(decoded.b_factor.as_ref(), atom, 0.0),
        formal_charge: (formal_charge, Presence::Present),
        atom_site_id: atom_site_id(decoded.atom_id.as_ref(), atom)?,
    })
}

/// Builds one normalized residue record from an MMTF group.
///
/// Converts component identifiers, sequence identifiers, insertion codes, and
/// entity classification into the normalized hierarchy representation.
fn residue_record(
    data: &mut StructureData,
    decoded: &Decoded,
    group: &Group,
    group_index: usize,
    chain_index: usize,
) -> Result<ResidueRecord, Diagnostic> {
    let component = intern(data, &group.name)?;

    let sequence = decoded
        .sequence_index
        .as_ref()
        .and_then(|values| values.get(group_index))
        .and_then(|value| value.checked_add(1))
        .map_or(OptionalI32::NONE, OptionalI32::some);

    let auth = decoded
        .group_id
        .get(group_index)
        .copied()
        .ok_or_else(|| schema_error("group id missing"))?;

    let ins_code = char_symbol(
        data,
        decoded
            .ins_code
            .as_ref()
            .and_then(|values| values.get(group_index))
            .copied(),
    )?
    .map_or(OptionalSymbol::NONE, OptionalSymbol::some);

    let kind = entity_kind_for_chain(decoded.file.entity_list.as_deref(), chain_index)?;

    Ok(ResidueRecord {
        label_comp_id: component,
        auth_comp_id: OptionalSymbol::some(component),
        label_seq_id: sequence,
        auth_seq_id: OptionalI32::some(auth),
        ins_code,
        het: kind != EntityKind::Polymer,
    })
}

/// Builds normalized entities and maps every MMTF chain to its owning entity.
///
/// Rejects missing entity metadata, invalid chain indices, duplicate ownership,
/// and chains that are not assigned to any entity.
fn build_entities(
    data: &mut StructureData,
    decoded: &Decoded,
) -> Result<Vec<EntityIndex>, Diagnostic> {
    let chain_count = usize_of(decoded.file.num_chains)?;
    let entities =
        decoded.file.entity_list.as_deref().ok_or_else(|| {
            schema_error("MMTF entityList is required by the normalized hierarchy")
        })?;

    let mut mapping = vec![None; chain_count];

    for (position, entity) in entities.iter().enumerate() {
        let ordinal = position
            .checked_add(1)
            .ok_or_else(|| schema_error("MMTF entity ordinal overflows"))?;
        let ordinal = ordinal.to_string();
        let id = intern(data, &ordinal)?;

        let description = if entity.description.is_empty() {
            OptionalSymbol::NONE
        } else {
            OptionalSymbol::some(intern(data, &entity.description)?)
        };

        let index = data
            .topology
            .entities
            .push(id, entity_kind(entity)?, description, &[])
            .map_err(|error| {
                schema_error("entity table rejected a row").with_context("cause", error.to_string())
            })?;

        for &chain in &entity.chain_index_list {
            let chain = usize_of(chain)?;
            let slot = mapping
                .get_mut(chain)
                .ok_or_else(|| schema_error("entity chain index is out of bounds"))?;

            if slot.replace(index).is_some() {
                return Err(schema_error("chain belongs to more than one MMTF entity"));
            }
        }
    }

    mapping
        .into_iter()
        .map(|entity| entity.ok_or_else(|| schema_error("entity mapping missing")))
        .collect()
}

/// Reads and validates an optional MMTF atom-site identifier.
///
/// An absent identifier list is represented by zero. When the source list is
/// present, identifiers must exist, be non-negative, and be non-zero.
fn atom_site_id(values: Option<&Vec<i32>>, atom: usize) -> Result<u32, Diagnostic> {
    let Some(values) = values else {
        return Ok(0);
    };

    let value = values
        .get(atom)
        .copied()
        .ok_or_else(|| schema_error("atom identifier missing"))?;

    let value = u32::try_from(value).map_err(|_| schema_error("atom identifier is negative"))?;

    if value == 0 {
        Err(schema_error(
            "present MMTF atom identifier cannot be represented as zero",
        ))
    } else {
        Ok(value)
    }
}

/// Adds bonds declared within a single MMTF group.
///
/// Bond endpoints are translated from group-local indices into retained
/// normalized atom indices. Bonds touching filtered atoms are ignored.
fn add_group_bonds(
    group: &Group,
    atom_offset: usize,
    keep_map: &[Option<AtomIndex>],
    bonds: &mut BondTableBuilder,
) -> Result<(), Diagnostic> {
    if !group.bond_atom_list.len().is_multiple_of(2) {
        return Err(schema_error("group bond endpoint list has odd length"));
    }

    for (bond_index, pair) in group.bond_atom_list.chunks_exact(2).enumerate() {
        let first = checked_atom_offset(atom_offset, pair[0])?;
        let second = checked_atom_offset(atom_offset, pair[1])?;

        let Some((atom_a, atom_b)) = retained_atoms(keep_map, first, second) else {
            continue;
        };

        bonds.push(BondRecord {
            atom_a,
            atom_b,
            order: decoded_order(
                &group.bond_order_list,
                &group.bond_resonance_list,
                bond_index,
            ),
            provenance: BondProvenance::File,
        });
    }

    Ok(())
}

/// Adds global MMTF bonds belonging to the selected model.
///
/// Global atom indices outside the model are ignored. Bonds touching atoms
/// that were filtered from the normalized structure are also ignored.
fn add_global_bonds(
    decoded: &Decoded,
    range: ModelRange,
    keep_map: &[Option<AtomIndex>],
    bonds: &mut BondTableBuilder,
) -> Result<(), Diagnostic> {
    if !decoded.global_bonds.len().is_multiple_of(2) {
        return Err(schema_error("global bond endpoint list has odd length"));
    }

    let atom_end = range
        .atom_start
        .checked_add(range.atom_count)
        .ok_or_else(|| schema_error("MMTF model atom range overflows"))?;

    for (bond_index, pair) in decoded.global_bonds.chunks_exact(2).enumerate() {
        let first = usize_of(pair[0])?;
        let second = usize_of(pair[1])?;

        if !index_in_range(first, range.atom_start, atom_end)
            || !index_in_range(second, range.atom_start, atom_end)
        {
            continue;
        }

        let Some((atom_a, atom_b)) = retained_atoms(keep_map, first, second) else {
            continue;
        };

        bonds.push(BondRecord {
            atom_a,
            atom_b,
            order: decoded_order(
                &decoded.global_orders,
                &decoded.global_resonance,
                bond_index,
            ),
            provenance: BondProvenance::File,
        });
    }

    Ok(())
}

/// Resolves two original atom indices into their retained normalized indices.
///
/// Returns `None` when either endpoint lies outside the mapping or corresponds
/// to an atom that was intentionally filtered out.
fn retained_atoms(
    keep_map: &[Option<AtomIndex>],
    first: usize,
    second: usize,
) -> Option<(AtomIndex, AtomIndex)> {
    match (keep_map.get(first), keep_map.get(second)) {
        (Some(Some(atom_a)), Some(Some(atom_b))) => Some((*atom_a, *atom_b)),
        _ => None,
    }
}

/// Converts a group-local signed atom index into an absolute atom index.
///
/// Rejects negative MMTF indices and arithmetic overflow while applying the
/// group atom offset.
fn checked_atom_offset(atom_offset: usize, local: i32) -> Result<usize, Diagnostic> {
    atom_offset
        .checked_add(usize_of(local)?)
        .ok_or_else(|| schema_error("MMTF atom index overflows"))
}

/// Returns whether an index belongs to a prevalidated half-open interval.
///
/// The caller supplies `end` after checked range construction, avoiding
/// repeated range arithmetic while processing bond endpoints.
fn index_in_range(index: usize, start: usize, end: usize) -> bool {
    index >= start && index < end
}

/// Decodes MMTF bond-order and resonance fields into the normalized bond order.
///
/// Resonance takes precedence and is represented as aromatic. Missing or
/// unsupported numeric bond orders are represented as unknown.
fn decoded_order(orders: &[i32], resonance: &[i32], index: usize) -> BondOrder {
    if resonance.get(index) == Some(&1) {
        return BondOrder::Aromatic;
    }

    match orders.get(index).copied() {
        Some(1) => BondOrder::Single,
        Some(2) => BondOrder::Double,
        Some(3) => BondOrder::Triple,
        Some(4) => BondOrder::Quadruple,
        _ => BondOrder::Unknown,
    }
}

/// Resolves the group template referenced by one decoded group instance.
///
/// Validates both the group-instance index and the referenced group-template
/// index before returning the template.
fn group_at(decoded: &Decoded, group_index: usize) -> Result<&Group, Diagnostic> {
    let type_index = decoded
        .group_type
        .get(group_index)
        .copied()
        .ok_or_else(|| schema_error("group type missing"))
        .and_then(usize_of)?;

    decoded
        .file
        .group_list
        .get(type_index)
        .ok_or_else(|| schema_error("group type index is out of bounds"))
}

/// Resolves the chemical element of a group-local atom.
///
/// Missing, out-of-range, or unrecognized element symbols map to
/// `Element::UNKNOWN`.
fn group_element(group: &Group, local: usize) -> Element {
    group
        .element_list
        .as_ref()
        .and_then(|elements| elements.get(local))
        .and_then(|symbol| Element::from_symbol(symbol))
        .map_or(Element::UNKNOWN, |element| element)
}

/// Counts atoms contained in a contiguous range of decoded MMTF groups.
///
/// Validates the requested group range, every referenced group template, and
/// all arithmetic used to accumulate the atom count.
fn atoms_in_groups(decoded: &Decoded, start: usize, count: usize) -> Result<usize, Diagnostic> {
    let end = start
        .checked_add(count)
        .ok_or_else(|| schema_error("MMTF group range overflows"))?;

    let mut atoms = 0usize;

    for group_index in start..end {
        let group_atoms = group_at(decoded, group_index)?.atom_name_list.len();
        atoms = atoms
            .checked_add(group_atoms)
            .ok_or_else(|| schema_error("MMTF atom count overflows"))?;
    }

    Ok(atoms)
}

include!("validation.rs");
