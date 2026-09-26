//! Direct encoding for canonical metadata categories.

use super::column::{FloatColumnBuilder, IntegerColumnBuilder, TextColumnBuilder};
use crate::container::EncodedCategory;
use molframe_cif::{CanonicalProjection, CanonicalValue};
use molframe_core::diagnostic::Diagnostic;
use molframe_core::topology::EntityKind;

pub(super) fn encode(
    projection: CanonicalProjection<'_>,
) -> Result<Vec<EncodedCategory>, Diagnostic> {
    let mut categories = Vec::new();
    encode_entry(projection, &mut categories)?;
    encode_cell(projection, &mut categories)?;
    encode_entities(projection, &mut categories)?;
    encode_references(projection, &mut categories)?;
    Ok(categories)
}

fn encode_entry(
    projection: CanonicalProjection<'_>,
    target: &mut Vec<EncodedCategory>,
) -> Result<(), Diagnostic> {
    let entry = &projection.structure().data().entry;
    for (category, item, value) in [
        ("entry", "id", entry.id.as_deref()),
        ("struct", "title", entry.title.as_deref()),
        ("exptl", "method", entry.method.as_deref()),
    ] {
        if let Some(value) = value {
            target.push(single_text(category, item, value)?);
        }
    }
    if let Some(value) = entry.resolution {
        let mut column = FloatColumnBuilder::new("ls_d_res_high", 1);
        column.push(CanonicalValue::Present(round(f64::from(value))));
        target.push(EncodedCategory {
            name: "_refine".to_owned(),
            row_count: 1,
            columns: vec![column.finish()?],
        });
    }
    Ok(())
}

fn encode_cell(
    projection: CanonicalProjection<'_>,
    target: &mut Vec<EncodedCategory>,
) -> Result<(), Diagnostic> {
    let Some(cell) = projection.structure().data().cell else {
        return Ok(());
    };
    let names = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ];
    let values = [
        cell.lengths[0],
        cell.lengths[1],
        cell.lengths[2],
        cell.angles[0],
        cell.angles[1],
        cell.angles[2],
    ];
    let mut columns = Vec::with_capacity(names.len());
    for (name, value) in names.into_iter().zip(values) {
        let mut column = FloatColumnBuilder::new(name, 1);
        column.push(CanonicalValue::Present(round(value)));
        columns.push(column.finish()?);
    }
    target.push(EncodedCategory {
        name: "_cell".to_owned(),
        row_count: 1,
        columns,
    });
    Ok(())
}

fn encode_entities(
    projection: CanonicalProjection<'_>,
    target: &mut Vec<EncodedCategory>,
) -> Result<(), Diagnostic> {
    let structure = projection.structure();
    let entities = &structure.data().topology.entities;
    if entities.is_empty() {
        return Ok(());
    }
    let mut ids = TextColumnBuilder::new("id", entities.len());
    let mut kinds = TextColumnBuilder::new("type", entities.len());
    let mut descriptions = TextColumnBuilder::new("pdbx_description", entities.len());
    for entity in entities.iter() {
        ids.push(inapplicable(
            entities
                .id(entity)
                .and_then(|value| structure.resolve(value)),
        ));
        kinds.push(entity_kind(entities.kind(entity)));
        descriptions.push(inapplicable(
            entities
                .description(entity)
                .and_then(|value| structure.resolve(value)),
        ));
    }
    target.push(EncodedCategory {
        name: "_entity".to_owned(),
        row_count: entities.len(),
        columns: vec![ids.finish()?, kinds.finish()?, descriptions.finish()?],
    });
    encode_entity_polymers(projection, target)?;
    encode_entity_sequence(projection, target)
}

/// Declares each polymer entity's type, so a reader recovers the chain's kind.
fn encode_entity_polymers(
    projection: CanonicalProjection<'_>,
    target: &mut Vec<EncodedCategory>,
) -> Result<(), Diagnostic> {
    let structure = projection.structure();
    let declared = molframe_cif::declared_polymer_types(structure);
    if declared.is_empty() {
        return Ok(());
    }
    let entities = &structure.data().topology.entities;
    let mut ids = TextColumnBuilder::new("entity_id", declared.len());
    let mut kinds = TextColumnBuilder::new("type", declared.len());
    for &(entity, kind) in &declared {
        ids.push(inapplicable(
            entities
                .id(entity)
                .and_then(|value| structure.resolve(value)),
        ));
        kinds.push(CanonicalValue::Present(kind));
    }
    target.push(EncodedCategory {
        name: "_entity_poly".to_owned(),
        row_count: declared.len(),
        columns: vec![ids.finish()?, kinds.finish()?],
    });
    Ok(())
}

fn encode_entity_sequence(
    projection: CanonicalProjection<'_>,
    target: &mut Vec<EncodedCategory>,
) -> Result<(), Diagnostic> {
    let structure = projection.structure();
    let entities = &structure.data().topology.entities;
    let rows = entities
        .iter()
        .map(|entity| entities.canonical_sequence(entity).len())
        .sum();
    if rows == 0 {
        return Ok(());
    }
    let mut ids = TextColumnBuilder::new("entity_id", rows);
    let mut numbers = IntegerColumnBuilder::new("num", rows);
    let mut components = TextColumnBuilder::new("mon_id", rows);
    for entity in entities.iter() {
        let id = entities
            .id(entity)
            .and_then(|value| structure.resolve(value));
        for (position, component) in entities.canonical_sequence(entity).iter().enumerate() {
            ids.push(inapplicable(id));
            let number = i64::try_from(position + 1)
                .map_or(CanonicalValue::Unknown, CanonicalValue::Present);
            numbers.push(number);
            components.push(inapplicable(structure.resolve(*component)));
        }
    }
    target.push(EncodedCategory {
        name: "_entity_poly_seq".to_owned(),
        row_count: rows,
        columns: vec![ids.finish()?, numbers.finish()?, components.finish()?],
    });
    Ok(())
}

fn encode_references(
    projection: CanonicalProjection<'_>,
    target: &mut Vec<EncodedCategory>,
) -> Result<(), Diagnostic> {
    let Some(references) = projection.sequence_references() else {
        return Ok(());
    };
    if !references.sequences.is_empty() {
        let rows = references.sequences.len();
        let mut columns: Vec<TextColumnBuilder> = [
            "id",
            "entity_id",
            "db_name",
            "db_code",
            "pdbx_db_accession",
            "pdbx_seq_one_letter_code",
        ]
        .map(|name| TextColumnBuilder::new(name, rows))
        .into();
        for value in &references.sequences {
            columns[0].push(CanonicalValue::Present(&value.id));
            columns[1].push(CanonicalValue::Present(&value.entity_id));
            columns[2].push(unknown(value.database_name.as_deref()));
            columns[3].push(unknown(value.database_code.as_deref()));
            columns[4].push(unknown(value.accession.as_deref()));
            columns[5].push(unknown(value.one_letter_code.as_deref()));
        }
        let columns = columns
            .into_iter()
            .map(TextColumnBuilder::finish)
            .collect::<Result<Vec<_>, _>>()?;
        target.push(EncodedCategory {
            name: "_struct_ref".to_owned(),
            row_count: rows,
            columns,
        });
    }
    encode_alignments(references, target)
}

fn encode_alignments(
    references: &molframe_core::structure::SequenceReferences,
    target: &mut Vec<EncodedCategory>,
) -> Result<(), Diagnostic> {
    if references.alignments.is_empty() {
        return Ok(());
    }
    let rows = references.alignments.len();
    let mut ids = TextColumnBuilder::new("align_id", rows);
    let mut reference_ids = TextColumnBuilder::new("ref_id", rows);
    let mut chains = TextColumnBuilder::new("pdbx_strand_id", rows);
    let mut numbers: Vec<IntegerColumnBuilder> = [
        "seq_align_beg",
        "seq_align_end",
        "db_align_beg",
        "db_align_end",
    ]
    .map(|name| IntegerColumnBuilder::new(name, rows))
    .into();
    for alignment in &references.alignments {
        ids.push(CanonicalValue::Present(&alignment.id));
        reference_ids.push(CanonicalValue::Present(&alignment.reference_id));
        let chain_ids = alignment.chain_ids.join(",");
        chains.push(unknown(
            (!chain_ids.is_empty()).then_some(chain_ids.as_str()),
        ));
        for (column, value) in numbers.iter_mut().zip([
            alignment.canonical[0],
            alignment.canonical[1],
            alignment.reference[0],
            alignment.reference[1],
        ]) {
            column.push(CanonicalValue::Present(i64::from(value)));
        }
    }
    let mut columns = vec![ids.finish()?, reference_ids.finish()?, chains.finish()?];
    columns.extend(
        numbers
            .into_iter()
            .map(IntegerColumnBuilder::finish)
            .collect::<Result<Vec<_>, _>>()?,
    );
    target.push(EncodedCategory {
        name: "_struct_ref_seq".to_owned(),
        row_count: rows,
        columns,
    });
    Ok(())
}

fn single_text(
    category: &str,
    item: &'static str,
    value: &str,
) -> Result<EncodedCategory, Diagnostic> {
    let mut column = TextColumnBuilder::new(item, 1);
    column.push(CanonicalValue::Present(value));
    Ok(EncodedCategory {
        name: format!("_{category}"),
        row_count: 1,
        columns: vec![column.finish()?],
    })
}

fn entity_kind(value: Option<EntityKind>) -> CanonicalValue<&'static str> {
    match value {
        Some(EntityKind::Polymer) => CanonicalValue::Present("polymer"),
        Some(EntityKind::NonPolymer) => CanonicalValue::Present("non-polymer"),
        Some(EntityKind::Water) => CanonicalValue::Present("water"),
        Some(EntityKind::Branched) => CanonicalValue::Present("branched"),
        Some(_) | None => CanonicalValue::Unknown,
    }
}

fn inapplicable(value: Option<&str>) -> CanonicalValue<&str> {
    match value.filter(|text| !text.is_empty()) {
        Some(text) => CanonicalValue::Present(text),
        None => CanonicalValue::Inapplicable,
    }
}

fn unknown(value: Option<&str>) -> CanonicalValue<&str> {
    match value.filter(|text| !text.is_empty()) {
        Some(text) => CanonicalValue::Present(text),
        None => CanonicalValue::Unknown,
    }
}

fn round(value: f64) -> f64 {
    (value * 10_000.0).round() / 10_000.0
}
