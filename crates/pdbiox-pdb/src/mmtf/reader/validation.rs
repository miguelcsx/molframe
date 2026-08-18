/// Validates that an MMTF version belongs to a supported major version.
///
/// Major versions zero and one are accepted; malformed or newer major versions
/// produce a schema diagnostic carrying the original version string.
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

/// Validates mandatory MMTF row counts against configured reader limits.
///
/// Every declared atom, group, chain, and model count must be positive,
/// representable as `usize`, and permitted by the configured row limit.
fn validate_declared_counts(file: &File, options: &ReadOptions) -> Result<(), Diagnostic> {
    for (name, value) in [
        ("numAtoms", file.num_atoms),
        ("numGroups", file.num_groups),
        ("numChains", file.num_chains),
        ("numModels", file.num_models),
    ] {
        validate_declared_count(name, value, options.limits.rows_per_category)?;
    }

    Ok(())
}

/// Validates one mandatory signed MMTF count.
///
/// Converts the count into a platform index, applies the configured row limit,
/// and rejects zero for fields required by the normalized hierarchy.
fn validate_declared_count(
    name: &'static str,
    value: i32,
    row_limit: u64,
) -> Result<(), Diagnostic> {
    let value = usize_of(value)?;

    if value as u64 > row_limit {
        return Err(pdbiox_core::io::Limits::exceeded("MMTF rows", value));
    }

    if value == 0 {
        return Err(schema_error("required MMTF count is zero").with_context("field", name));
    }

    Ok(())
}

/// Verifies decoded array lengths against the counts declared by the MMTF file.
///
/// Coordinates, group fields, chain identifiers, and group-per-chain entries
/// must exactly match their corresponding declared dimensions.
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
        validate_exact_length(name, actual, expected)?;
    }

    Ok(())
}

/// Verifies that one decoded MMTF field has its required number of values.
///
/// Length mismatches include the field name and both dimensions in the emitted
/// schema diagnostic.
fn validate_exact_length(
    name: &'static str,
    actual: usize,
    expected: usize,
) -> Result<(), Diagnostic> {
    if actual == expected {
        return Ok(());
    }

    Err(
        schema_error("MMTF field length disagrees with declared count")
            .with_context("field", name)
            .with_context("expected", expected.to_string())
            .with_context("actual", actual.to_string()),
    )
}

/// Reads one Cartesian coordinate triple from the decoded MMTF arrays.
///
/// Each axis is independently bounds-checked so truncated coordinate arrays
/// yield a schema diagnostic instead of a panic.
fn position(decoded: &Decoded, atom: usize) -> Result<[f32; 3], Diagnostic> {
    let x = decoded
        .x
        .get(atom)
        .copied()
        .ok_or_else(|| schema_error("x coordinate missing"))?;
    let y = decoded
        .y
        .get(atom)
        .copied()
        .ok_or_else(|| schema_error("y coordinate missing"))?;
    let z = decoded
        .z
        .get(atom)
        .copied()
        .ok_or_else(|| schema_error("z coordinate missing"))?;

    Ok([x, y, z])
}

/// Resolves an optional per-atom floating-point field.
///
/// Missing lists or missing entries yield the supplied default and mark the
/// field as inapplicable; existing entries retain their source value.
fn optional_float(values: Option<&Vec<f32>>, atom: usize, default: f32) -> (f32, Presence) {
    match values.and_then(|values| values.get(atom)).copied() {
        Some(value) => (value, Presence::Present),
        None => (default, Presence::Inapplicable),
    }
}

/// Converts an optional six-value MMTF unit cell into normalized precision.
///
/// The source sequence must contain exactly three lengths followed by three
/// angles when present.
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

/// Converts an MMTF entity type into the normalized entity classification.
///
/// Known MMTF aliases are normalized, while unsupported entity types are
/// reported with their original textual value.
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

/// Determines the entity classification assigned to an MMTF chain.
///
/// A missing entity list or a chain not referenced by any entity is represented
/// as `EntityKind::Unknown`.
fn entity_kind_for_chain(
    entities: Option<&[Entity]>,
    chain: usize,
) -> Result<EntityKind, Diagnostic> {
    let Some(entities) = entities else {
        return Ok(EntityKind::Unknown);
    };

    let Ok(chain) = i32::try_from(chain) else {
        return Ok(EntityKind::Unknown);
    };

    entities
        .iter()
        .find(|entity| entity.chain_index_list.contains(&chain))
        .map_or(Ok(EntityKind::Unknown), entity_kind)
}

/// Converts an optional encoded MMTF character into an interned symbol.
///
/// Zero and absent values represent no symbol. ASCII values are interned
/// directly from a stack buffer; non-ASCII bytes are rejected.
fn char_symbol(
    data: &mut StructureData,
    value: Option<u8>,
) -> Result<Option<SymbolId>, Diagnostic> {
    match value {
        Some(0) | None => Ok(None),
        Some(value) if value.is_ascii() => {
            let mut buffer = [0u8; 4];
            let text = char::from(value).encode_utf8(&mut buffer);
            intern(data, text).map(Some)
        }
        Some(_) => Err(schema_error("MMTF character is not ASCII")),
    }
}

/// Interns an identifier into the structure dictionary.
///
/// Dictionary-capacity failures are translated into the reader's standard
/// resource-limit diagnostic.
fn intern(data: &mut StructureData, text: &str) -> Result<SymbolId, Diagnostic> {
    data.dictionary
        .intern(text)
        .map_err(|_| pdbiox_core::io::Limits::exceeded("identifier dictionary", text))
}

/// Converts a non-negative MMTF count or index into `usize`.
///
/// Negative source values are treated as malformed schema data rather than
/// being reinterpreted through a lossy cast.
fn usize_of(value: i32) -> Result<usize, Diagnostic> {
    usize::try_from(value).map_err(|_| schema_error("MMTF count or index is negative"))
}

/// Creates the standard schema diagnostic used for malformed MMTF input.
///
/// The supplied static message is attached to the parser's canonical schema
/// error code.
fn schema_error(message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1102).with_message(message)
}

