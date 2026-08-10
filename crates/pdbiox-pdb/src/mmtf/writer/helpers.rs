fn encoded_order(order: BondOrder) -> (i32, i32) {
    match order {
        BondOrder::Single | BondOrder::Polymeric => (1, 0),
        BondOrder::Double => (2, 0),
        BondOrder::Triple => (3, 0),
        BondOrder::Quadruple => (4, 0),
        BondOrder::Aromatic => (1, 1),
        BondOrder::Unknown => (-1, -1),
    }
}

type OptionalEncodedBonds = (
    Option<serde_bytes::ByteBuf>,
    Option<serde_bytes::ByteBuf>,
    Option<serde_bytes::ByteBuf>,
);

fn encode_bond_columns(
    atoms: &[i32],
    orders: &[i32],
    resonance: &[i32],
) -> Result<OptionalEncodedBonds, Diagnostic> {
    Ok((
        (!atoms.is_empty()).then(|| encode_i32(atoms)).transpose()?,
        (!orders.is_empty())
            .then(|| encode_i32(orders))
            .transpose()?,
        (!resonance.is_empty())
            .then(|| encode_i32(resonance))
            .transpose()?,
    ))
}

fn required_text(value: Option<&str>, field: &'static str) -> Result<String, Diagnostic> {
    let value = value.ok_or_else(|| refusal("MMTF required identifier is absent"))?;
    if !value.is_ascii() {
        return Err(refusal("MMTF identifiers must be ASCII").with_context("field", field));
    }
    Ok(value.to_string())
}

fn single_ascii(value: Option<&str>, field: &'static str) -> Result<u8, Diagnostic> {
    let Some(value) = value else { return Ok(0) };
    let [byte] = value.as_bytes() else {
        return Err(refusal("MMTF character field must contain one ASCII byte")
            .with_context("field", field));
    };
    if !byte.is_ascii() {
        return Err(refusal("MMTF character field must be ASCII").with_context("field", field));
    }
    Ok(*byte)
}

fn i32_of(value: usize, field: &'static str) -> Result<i32, Diagnostic> {
    i32::try_from(value)
        .map_err(|_| refusal("MMTF count exceeds signed 32-bit range").with_context("field", field))
}

fn source_metadata(structure: &Structure) -> Result<&MmtfMetadata, Diagnostic> {
    structure
        .extensions()
        .get::<MmtfMetadata>(MMTF_METADATA_EXTENSION)
        .ok_or_else(|| refusal("MMTF chemistry metadata is absent"))
}

fn encoded_cell(structure: &Structure) -> Result<Option<Vec<f32>>, Diagnostic> {
    structure.data().cell.map_or(Ok(None), |cell| {
        cell.lengths
            .into_iter()
            .chain(cell.angles)
            .map(|value| {
                value.to_f32().ok_or_else(|| {
                    refusal("unit-cell value exceeds the MMTF floating-point range")
                })
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some)
    })
}

fn experimental_methods(structure: &Structure) -> Option<Vec<String>> {
    structure
        .data()
        .entry
        .method
        .as_deref()
        .map(|method| vec![method.to_string()])
}

fn required_presence<'a, T>(
    values: Option<&'a [T]>,
    expected: bool,
    field: &'static str,
) -> Result<Option<&'a [T]>, Diagnostic> {
    match (values, expected) {
        (Some(values), true) => Ok(Some(values)),
        (None, false) => Ok(None),
        _ => Err(refusal("MMTF optional field presence changed").with_context("field", field)),
    }
}

fn encode_optional_f32(
    values: Option<&[f32]>,
    expected: bool,
    field: &'static str,
) -> Result<Option<serde_bytes::ByteBuf>, Diagnostic> {
    required_presence(values, expected, field)?
        .map(encode_f32)
        .transpose()
}

fn encode_optional_i32(
    values: Option<&[i32]>,
    expected: bool,
    field: &'static str,
) -> Result<Option<serde_bytes::ByteBuf>, Diagnostic> {
    required_presence(values, expected, field)?
        .map(encode_i32)
        .transpose()
}

fn chars_with_presence(
    values: &[u8],
    expected: bool,
    field: &'static str,
) -> Result<Option<serde_bytes::ByteBuf>, Diagnostic> {
    if expected {
        Ok(Some(encode_chars(values)?))
    } else if values.iter().all(|value| *value == 0) {
        Ok(None)
    } else {
        Err(refusal("MMTF optional field presence changed").with_context("field", field))
    }
}

fn validate_elements(
    atoms: &[pdbiox_core::structure::AtomRef<'_>],
    elements: &[Box<str>],
) -> Result<Vec<String>, Diagnostic> {
    if atoms.len() != elements.len()
        || atoms.iter().zip(elements).any(|(atom, deposited)| {
            atom.element()
                .is_some_and(|element| element.symbol() != deposited.as_ref())
        })
    {
        return Err(refusal(
            "MMTF deposited elements no longer match the structure",
        ));
    }
    Ok(elements.iter().map(ToString::to_string).collect())
}

fn refusal(message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4105).with_message(message)
}
