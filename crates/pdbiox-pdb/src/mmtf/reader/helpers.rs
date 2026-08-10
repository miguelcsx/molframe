fn same_topology(left: &Structure, right: &Structure) -> bool {
    if left.atom_count() != right.atom_count()
        || left.residue_count() != right.residue_count()
        || left.chain_count() != right.chain_count()
        || left.data().bonds.iter().ne(right.data().bonds.iter())
    {
        return false;
    }
    let chains_match = left
        .data()
        .chains()
        .zip(right.data().chains())
        .all(|(a, b)| {
            a.label() == b.label()
                && a.auth_label() == b.auth_label()
                && a.polymer_kind() == b.polymer_kind()
        });
    let residues_match = left
        .data()
        .residues()
        .zip(right.data().residues())
        .all(|(a, b)| {
            a.name() == b.name()
                && a.auth_seq_id() == b.auth_seq_id()
                && a.label_seq_id() == b.label_seq_id()
                && a.ins_code() == b.ins_code()
                && a.is_het() == b.is_het()
        });
    let atoms_match = left.data().atoms().zip(right.data().atoms()).all(|(a, b)| {
        a.name() == b.name()
            && a.element() == b.element()
            && a.alt_label() == b.alt_label()
            && a.formal_charge() == b.formal_charge()
    });
    chains_match && residues_match && atoms_match
}

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
    let charge = *group
        .formal_charge_list
        .get(local)
        .ok_or_else(|| schema_error("group formal charge missing"))?;
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
        formal_charge: (
            i8::try_from(charge).map_err(|_| schema_error("formal charge is outside i8"))?,
            Presence::Present,
        ),
        atom_site_id: atom_site_id(decoded.atom_id.as_ref(), atom)?,
    })
}

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
    let auth = *decoded
        .group_id
        .get(group_index)
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

fn build_entities(
    data: &mut StructureData,
    decoded: &Decoded,
) -> Result<Vec<EntityIndex>, Diagnostic> {
    let chain_count = usize_of(decoded.file.num_chains)?;
    let mut mapping = vec![None; chain_count];
    let entities =
        decoded.file.entity_list.as_deref().ok_or_else(|| {
            schema_error("MMTF entityList is required by the normalized hierarchy")
        })?;
    for (position, entity) in entities.iter().enumerate() {
        let id = intern(data, &(position + 1).to_string())?;
        let description = (!entity.description.is_empty())
            .then(|| intern(data, &entity.description))
            .transpose()?
            .map_or(OptionalSymbol::NONE, OptionalSymbol::some);
        let index = data
            .topology
            .entities
            .push(id, entity_kind(entity)?, description, &[])
            .map_err(|error| schema_error("entity table rejected a row").with_context("cause", error.to_string()))?;
        for chain in &entity.chain_index_list {
            let chain = usize_of(*chain)?;
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
        .map(|value| value.ok_or_else(|| schema_error("entity mapping missing")))
        .collect()
}

fn atom_site_id(values: Option<&Vec<i32>>, atom: usize) -> Result<u32, Diagnostic> {
    let Some(values) = values else {
        return Ok(0);
    };
    let value = *values
        .get(atom)
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
        let first = atom_offset + usize_of(pair[0])?;
        let second = atom_offset + usize_of(pair[1])?;
        let (Some(Some(atom_a)), Some(Some(atom_b))) = (keep_map.get(first), keep_map.get(second))
        else {
            continue;
        };
        bonds.push(BondRecord {
            atom_a: *atom_a,
            atom_b: *atom_b,
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

fn add_global_bonds(
    decoded: &Decoded,
    range: ModelRange,
    keep_map: &[Option<AtomIndex>],
    bonds: &mut BondTableBuilder,
) -> Result<(), Diagnostic> {
    if !decoded.global_bonds.len().is_multiple_of(2) {
        return Err(schema_error("global bond endpoint list has odd length"));
    }
    for (bond_index, pair) in decoded.global_bonds.chunks_exact(2).enumerate() {
        let (first, second) = (usize_of(pair[0])?, usize_of(pair[1])?);
        let atoms = range.atom_start..range.atom_start + range.atom_count;
        if !atoms.contains(&first) || !atoms.contains(&second) {
            continue;
        }
        let (Some(atom_a), Some(atom_b)) = (keep_map[first], keep_map[second]) else {
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

fn decoded_order(orders: &[i32], resonance: &[i32], index: usize) -> BondOrder {
    if resonance.get(index) == Some(&1) {
        return BondOrder::Aromatic;
    }
    let order = match orders.get(index) {
        Some(order) => *order,
        None => -1,
    };
    match order {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Quadruple,
        _ => BondOrder::Unknown,
    }
}

fn group_at(decoded: &Decoded, group_index: usize) -> Result<&Group, Diagnostic> {
    let type_index = usize_of(
        *decoded
            .group_type
            .get(group_index)
            .ok_or_else(|| schema_error("group type missing"))?,
    )?;
    decoded
        .file
        .group_list
        .get(type_index)
        .ok_or_else(|| schema_error("group type index is out of bounds"))
}

fn group_element(group: &Group, local: usize) -> Element {
    match group
        .element_list
        .as_ref()
        .and_then(|elements| elements.get(local))
        .and_then(|symbol| Element::from_symbol(symbol))
    {
        Some(element) => element,
        None => Element::UNKNOWN,
    }
}

fn atoms_in_groups(decoded: &Decoded, start: usize, count: usize) -> Result<usize, Diagnostic> {
    let end = start
        .checked_add(count)
        .ok_or_else(|| schema_error("MMTF group range overflows"))?;
    (start..end).try_fold(0usize, |sum, group| {
        sum.checked_add(group_at(decoded, group)?.atom_name_list.len())
            .ok_or_else(|| schema_error("MMTF atom count overflows"))
    })
}

fn validate_version(version: &str) -> Result<(), Diagnostic> {
    let major = version
        .split('.')
        .next()
        .and_then(|value| value.parse::<u32>().ok());
    if matches!(major, Some(0 | 1)) {
        Ok(())
    } else {
        Err(schema_error("unsupported MMTF major version").with_context("version", version))
    }
}

fn validate_declared_counts(file: &File, options: &ReadOptions) -> Result<(), Diagnostic> {
    for (name, value) in [
        ("numAtoms", file.num_atoms),
        ("numGroups", file.num_groups),
        ("numChains", file.num_chains),
        ("numModels", file.num_models),
    ] {
        let value = usize_of(value)?;
        if value as u64 > options.limits.rows_per_category {
            return Err(pdbiox_core::io::Limits::exceeded("MMTF rows", value));
        }
        if value == 0 {
            return Err(schema_error("required MMTF count is zero").with_context("field", name));
        }
    }
    Ok(())
}

fn validate_decoded_lengths(decoded: &Decoded) -> Result<(), Diagnostic> {
    let atoms = usize_of(decoded.file.num_atoms)?;
    let groups = usize_of(decoded.file.num_groups)?;
    let chains = usize_of(decoded.file.num_chains)?;
    for (name, actual, expected) in [
        ("xCoordList", decoded.x.len(), atoms),
        ("yCoordList", decoded.y.len(), atoms),
        ("zCoordList", decoded.z.len(), atoms),
        ("groupIdList", decoded.group_id.len(), groups),
        ("groupTypeList", decoded.group_type.len(), groups),
        ("chainIdList", decoded.chain_id.len(), chains),
        (
            "groupsPerChain",
            decoded.file.groups_per_chain.len(),
            chains,
        ),
    ] {
        if actual != expected {
            return Err(
                schema_error("MMTF field length disagrees with declared count")
                    .with_context("field", name)
                    .with_context("expected", expected.to_string())
                    .with_context("actual", actual.to_string()),
            );
        }
    }
    Ok(())
}

fn position(decoded: &Decoded, atom: usize) -> Result<[f32; 3], Diagnostic> {
    Ok([
        *decoded
            .x
            .get(atom)
            .ok_or_else(|| schema_error("x coordinate missing"))?,
        *decoded
            .y
            .get(atom)
            .ok_or_else(|| schema_error("y coordinate missing"))?,
        *decoded
            .z
            .get(atom)
            .ok_or_else(|| schema_error("z coordinate missing"))?,
    ])
}

fn optional_float(values: Option<&Vec<f32>>, atom: usize, default: f32) -> (f32, Presence) {
    values
        .and_then(|values| values.get(atom).copied())
        .map_or((default, Presence::Inapplicable), |value| {
            (value, Presence::Present)
        })
}

fn unit_cell(values: Option<&[f32]>) -> Result<Option<UnitCell>, Diagnostic> {
    let Some(values) = values else {
        return Ok(None);
    };
    let [a, b, c, alpha, beta, gamma] = values else {
        return Err(schema_error("unitCell must contain six values"));
    };
    Ok(Some(UnitCell {
        lengths: [f64::from(*a), f64::from(*b), f64::from(*c)],
        angles: [f64::from(*alpha), f64::from(*beta), f64::from(*gamma)],
    }))
}

fn entity_kind(entity: &Entity) -> Result<EntityKind, Diagnostic> {
    match entity.kind.as_str() {
        "polymer" => Ok(EntityKind::Polymer),
        "non-polymer" | "macrolide" => Ok(EntityKind::NonPolymer),
        "water" => Ok(EntityKind::Water),
        "branched" => Ok(EntityKind::Branched),
        value => Err(schema_error("unsupported MMTF entity type")
            .with_context("entity_type", value.to_owned())),
    }
}

fn entity_kind_for_chain(
    entities: Option<&[Entity]>,
    chain: usize,
) -> Result<EntityKind, Diagnostic> {
    let entities = match entities {
        Some(entities) => entities,
        None => &[],
    };
    entities
        .iter()
        .find(|entity| {
            entity
                .chain_index_list
                .iter()
                .any(|index| usize::try_from(*index) == Ok(chain))
        })
        .map_or(Ok(EntityKind::Unknown), entity_kind)
}

fn char_symbol(
    data: &mut StructureData,
    value: Option<u8>,
) -> Result<Option<SymbolId>, Diagnostic> {
    match value {
        Some(0) | None => Ok(None),
        Some(value) if value.is_ascii() => intern(data, &char::from(value).to_string()).map(Some),
        Some(_) => Err(schema_error("MMTF character is not ASCII")),
    }
}

fn intern(data: &mut StructureData, text: &str) -> Result<SymbolId, Diagnostic> {
    data.dictionary
        .intern(text)
        .map_err(|_| pdbiox_core::io::Limits::exceeded("identifier dictionary", text))
}

fn usize_of(value: i32) -> Result<usize, Diagnostic> {
    usize::try_from(value).map_err(|_| schema_error("MMTF count or index is negative"))
}

fn schema_error(message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1102).with_message(message)
}
