/// Moves the single coordinate frame out of a just-built model.
fn take_frame(data: &mut StructureData) -> Result<CoordinateBlock, Diagnostic> {
    let coordinates = std::mem::replace(
        &mut data.coords,
        CoordinateStore::Single(CoordinateBlock::new()),
    );
    match coordinates {
        CoordinateStore::Single(frame) => Ok(frame),
        CoordinateStore::Dense { .. } | CoordinateStore::Ragged { .. } => {
            Err(schema_error("MMTF model builder produced an ensemble"))
        }
    }
}

/// Finishes a dense ensemble without copying any coordinate frame.
fn dense_ensemble(
    mut template: StructureData,
    frames: Vec<CoordinateBlock>,
) -> Result<Structure, Diagnostic> {
    template.topology.models = pdbiox_core::topology::ModelTable::default();
    let chain_count = u32_of(template.topology.chains.len(), "chain count")?;
    for model in 0..frames.len() {
        push_model(&mut template, model, chain_count)?;
    }
    template.coords = CoordinateStore::Dense { frames };
    Ok(Structure::new(template))
}

/// Switches from the tentative dense path to ragged storage at the first
/// topology mismatch. Earlier equal models reuse the template's copy-on-write
/// topology and move their coordinate frames, so no model is reparsed.
fn build_ragged_after_mismatch(
    context: DecodeContext<'_>,
    model_bonds: &[Vec<GlobalBond>],
    ranges: &[ModelRange],
    template: &StructureData,
    frames: Vec<CoordinateBlock>,
    mismatch: usize,
    candidate: StructureData,
) -> Result<Structure, Diagnostic> {
    let mut models = Vec::with_capacity(ranges.len());
    for (model, frame) in frames.into_iter().enumerate() {
        models.push(model_from_template(template, frame, model)?);
    }
    models.push(Structure::new(candidate));
    for (model, &range) in ranges.iter().enumerate().skip(mismatch + 1) {
        models.push(Structure::new(build_model_data(
            context.decoded,
            context.options,
            context.base,
            context.entity_by_chain,
            &model_bonds[model],
            model,
            range,
        )?));
    }
    ragged_ensemble(models)
}

/// Creates one independent ragged model over a shared copy-on-write topology.
fn model_from_template(
    template: &StructureData,
    frame: CoordinateBlock,
    model: usize,
) -> Result<Structure, Diagnostic> {
    let mut data = template.clone();
    data.topology.models = pdbiox_core::topology::ModelTable::default();
    let chain_count = u32_of(data.topology.chains.len(), "chain count")?;
    push_model(&mut data, model, chain_count)?;
    data.coords = CoordinateStore::Single(frame);
    Ok(Structure::new(data))
}

/// Adds one model row with validated indices.
fn push_model(
    data: &mut StructureData,
    model: usize,
    chain_count: u32,
) -> Result<(), Diagnostic> {
    data.topology
        .models
        .push(
            i32::try_from(model + 1).map_err(|_| schema_error("too many models"))?,
            0..chain_count,
        )
        .map(|_| ())
        .map_err(|error| {
            schema_error("model table rejected a row").with_context("cause", error.to_string())
        })
}
