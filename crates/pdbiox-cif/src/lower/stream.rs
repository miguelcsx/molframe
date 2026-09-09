//! Final assembly for coordinate models lowered while CIF text is parsed.

use super::atoms::{AtomBuilder, AtomSiteRow, AtomSiteRowSink};
use super::entry::{assemble_ragged, finish_model, prepare_model};
use super::metadata::refresh;
use super::ragged::RaggedParts;
use crate::document::Document;
use pdbiox_core::coords::CoordinateBlock;
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::io::{ReadOptions, ReadResult};
use pdbiox_core::structure::{CoordinateStore, StructureData};
use pdbiox_core::topology::ModelTable;

pub(crate) struct StreamModelParts {
    pub(crate) number: i64,
    pub(crate) data: StructureData,
    pub(crate) findings: Diagnostics,
    pub(crate) coordinates: CoordinateBlock,
}

pub(crate) struct StreamFrameParts {
    pub(crate) number: i64,
    pub(crate) findings: Diagnostics,
    pub(crate) coordinates: CoordinateBlock,
}

pub(crate) enum StreamedModels {
    Dense {
        first: Box<StreamModelParts>,
        additional: Vec<StreamFrameParts>,
    },
    Ragged(Vec<StreamModelParts>),
}

/// Builds one model at a time, comparing later rows to the first model in place.
pub(crate) struct StreamModelBuilder<'a> {
    template: StructureData,
    asym_entities: Vec<super::entry::AsymEntity>,
    options: &'a ReadOptions,
    current: Option<AtomBuilder<'a>>,
    current_number: Option<i64>,
    first_capacity: usize,
    reserve_first: bool,
}

impl<'a> StreamModelBuilder<'a> {
    pub(crate) fn new(
        metadata: &Document,
        options: &'a ReadOptions,
        first_capacity: usize,
    ) -> Result<Self, Vec<Diagnostic>> {
        let Some(block) = metadata.first_block() else {
            return Err(vec![Diagnostic::new(Code::E1106)]);
        };
        let (template, _findings, asym_entities) = prepare_model(block, options, Vec::new());
        Ok(Self {
            template,
            asym_entities,
            options,
            current: None,
            current_number: None,
            first_capacity,
            reserve_first: true,
        })
    }

    pub(crate) fn start(&mut self, number: i64, reference: Option<&StructureData>) {
        let builder = AtomBuilder::new(
            self.template.clone(),
            Diagnostics::new(),
            self.options,
            self.asym_entities.clone(),
        );
        let mut builder = match reference {
            Some(reference) => builder.against(reference),
            None => builder.without_identity_tracking(),
        };
        if self.reserve_first {
            builder.reserve_atoms(self.first_capacity);
            self.reserve_first = false;
        } else if let Some(reference) = reference {
            builder.reserve_atoms(reference.topology.atom_count() as usize);
        }
        self.current = Some(builder);
        self.current_number = Some(number);
    }

    pub(crate) fn feed(&mut self, row: &dyn AtomSiteRow) {
        if let Some(builder) = &mut self.current {
            builder.feed(row);
        }
    }

    pub(crate) fn close(&mut self) -> Option<(StreamModelParts, bool)> {
        let builder = self.current.take()?;
        let number = self.current_number.take()?;
        let matches_reference = builder.matches_expected();
        let (data, findings, coordinates) = builder.finish();
        let coordinates = match coordinates {
            CoordinateStore::Single(block) => block,
            CoordinateStore::Dense { mut frames } => match frames.pop() {
                Some(block) => block,
                None => CoordinateBlock::new(),
            },
            CoordinateStore::Ragged { .. } => CoordinateBlock::new(),
        };
        Some((
            StreamModelParts {
                number,
                data,
                findings,
                coordinates,
            },
            matches_reference,
        ))
    }
}

pub(crate) fn share_model_topology(
    first: &StreamModelParts,
    frame: StreamFrameParts,
) -> StreamModelParts {
    let mut data = first.data.clone();
    let mut models = ModelTable::default();
    let chains = match u32::try_from(data.topology.chains.len()) {
        Ok(count) => 0..count,
        Err(_) => 0..0,
    };
    let number = match i32::try_from(frame.number) {
        Ok(number) => number,
        Err(_) => i32::MAX,
    };
    let _model_was_not_representable = models.push(number, chains).is_err();
    data.topology.models = models;
    StreamModelParts {
        number: frame.number,
        data,
        findings: frame.findings,
        coordinates: frame.coordinates,
    }
}

pub(crate) fn finish_streamed(
    metadata: &Document,
    options: &ReadOptions,
    parser_findings: Vec<Diagnostic>,
    models: StreamedModels,
) -> ReadResult {
    let Some(block) = metadata.first_block() else {
        return Err(vec![Diagnostic::new(Code::E1106)]);
    };
    match models {
        StreamedModels::Dense { first, additional } => {
            finish_dense(block, options, parser_findings, *first, additional)
        }
        StreamedModels::Ragged(models) => finish_ragged(block, options, parser_findings, models),
    }
}

fn finish_dense(
    block: &crate::document::DataBlock,
    options: &ReadOptions,
    parser_findings: Vec<Diagnostic>,
    first: StreamModelParts,
    additional: Vec<StreamFrameParts>,
) -> ReadResult {
    let mut findings = Diagnostics::with_capacity(parser_findings.len() + first.findings.len());
    findings.extend(parser_findings);
    findings.extend(first.findings.finish());

    let mut data = first.data;
    refresh(block, options, &mut data, &mut findings);
    let mut numbers = Vec::with_capacity(additional.len() + 1);
    let mut frames = Vec::with_capacity(additional.len() + 1);
    numbers.push(first.number);
    frames.push(first.coordinates);
    for frame in additional {
        numbers.push(frame.number);
        frames.push(frame.coordinates);
        findings.extend(frame.findings.finish());
    }
    replace_dense_models(&mut data, &numbers, &mut findings);
    let coordinates = if frames.len() == 1 {
        CoordinateStore::Single(match frames.pop() {
            Some(frame) => frame,
            None => CoordinateBlock::new(),
        })
    } else {
        CoordinateStore::Dense { frames }
    };
    let (structure, findings) = finish_model(block, data, findings, coordinates);
    options.finish(structure, findings)
}

fn replace_dense_models(data: &mut StructureData, numbers: &[i64], findings: &mut Diagnostics) {
    let chain_count = if let Ok(count) = u32::try_from(data.topology.chains.len()) {
        count
    } else {
        findings.push(Diagnostic::new(Code::E1901));
        0
    };
    let mut models = ModelTable::default();
    for number in numbers {
        let deposited = if let Ok(number) = i32::try_from(*number) {
            number
        } else {
            findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("model number is outside the supported integer range")
                    .in_category("atom_site"),
            );
            i32::MAX
        };
        if models.push(deposited, 0..chain_count).is_err() {
            findings.push(Diagnostic::new(Code::E3001));
        }
    }
    data.topology.models = models;
}

fn finish_ragged(
    block: &crate::document::DataBlock,
    options: &ReadOptions,
    parser_findings: Vec<Diagnostic>,
    models: Vec<StreamModelParts>,
) -> ReadResult {
    let mut findings = Diagnostics::with_capacity(parser_findings.len());
    findings.extend(parser_findings);
    let mut structures = Vec::with_capacity(models.len());
    let mut numbers = Vec::with_capacity(models.len());
    for model in models {
        let mut model_findings = Diagnostics::with_capacity(model.findings.len());
        model_findings.extend(model.findings.finish());
        let mut data = model.data;
        refresh(block, options, &mut data, &mut model_findings);
        let (structure, model_findings) = finish_model(
            block,
            data,
            model_findings,
            CoordinateStore::Single(model.coordinates),
        );
        findings.extend(model_findings);
        structures.push(structure);
        numbers.push(model.number);
    }
    let (structure, findings) = assemble_ragged(
        block,
        options,
        RaggedParts {
            models: structures,
            model_numbers: numbers,
            findings,
        },
    );
    options.finish(structure, findings)
}
