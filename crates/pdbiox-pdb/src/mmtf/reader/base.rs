/// Builds metadata shared by every model and interns the complete finite set
/// of identifiers once. Model builders then clone only copy-on-write handles;
/// repeated symbol lookups do not detach the shared dictionary.
fn base_data(decoded: &Decoded) -> Result<(StructureData, Vec<EntityIndex>), Diagnostic> {
    let mut data = StructureData::empty();
    data.entry.id = decoded.file.structure_id.clone().map(Into::into);
    data.entry.title = decoded.file.title.clone().map(Into::into);
    data.entry.method = decoded
        .file
        .experimental_methods
        .as_ref()
        .map(|methods| methods.join(", ").into_boxed_str());
    data.entry.resolution = decoded.file.resolution;
    data.cell = unit_cell(decoded.file.unit_cell.as_deref())?;
    let entity_by_chain = build_entities(&mut data, decoded)?;
    intern_model_identifiers(&mut data, decoded)?;
    Ok((data, entity_by_chain))
}

/// Interns identifiers referenced by any model before the dictionary is shared.
fn intern_model_identifiers(
    data: &mut StructureData,
    decoded: &Decoded,
) -> Result<(), Diagnostic> {
    for group in &decoded.file.group_list {
        intern(data, &group.name)?;
        for atom_name in &group.atom_name_list {
            intern(data, atom_name)?;
        }
    }
    for chain in &decoded.chain_id {
        intern(data, chain)?;
    }
    if let Some(chains) = &decoded.chain_name {
        for chain in chains {
            intern(data, chain)?;
        }
    }

    let mut characters = [false; 128];
    for value in decoded
        .alt_loc
        .iter()
        .chain(&decoded.ins_code)
        .flat_map(|values| values.iter().copied())
        .filter(|value| *value != 0 && value.is_ascii())
    {
        characters[usize::from(value)] = true;
    }
    for (value, present) in characters.into_iter().enumerate() {
        if present {
            let value = u8::try_from(value)
                .map_err(|_| schema_error("MMTF character index exceeds byte range"))?;
            char_symbol(data, Some(value))?;
        }
    }
    Ok(())
}
