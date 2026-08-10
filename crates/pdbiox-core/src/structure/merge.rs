//! Deterministic concatenation of compatible immutable structures.

use super::merge_annotations::merge_annotations;
use super::{CoordinateStore, EntryMetadata, Structure, StructureData, UnitCell, validate};
use crate::bond::{BondRecord, BondTableBuilder};
use crate::chunk::{AtomRecord, ChunkBuilder};
use crate::coords::{CoordinateBlock, CoordinateGeneration};
use crate::diagnostic::{Code, Diagnostic};
use crate::index::{AtomIndex, EntityIndex, ModelIndex, ResidueIndex};
use crate::optional::OptionalSymbol;
use crate::symbol::{AltId, Interner, SymbolId};
use crate::topology::{ChainRecord, ResidueRecord, Topology};

impl Structure {
    /// Concatenates compatible structures into one immutable system.
    ///
    /// Atom, residue, chain and entity order follows input order. Dense inputs
    /// must have equal frame counts; their corresponding frames are
    /// concatenated. Dictionary identifiers, bonds and symbolic annotations are
    /// remapped rather than compared by structure-local integer value.
    ///
    /// # Errors
    ///
    /// Returns diagnostics for invalid inputs, ragged or unequal model axes,
    /// incompatible same-named annotation types, or representational overflow.
    pub fn merge(structures: &[Self]) -> Result<Self, Vec<Diagnostic>> {
        if structures.is_empty() {
            return Ok(Self::new(StructureData::empty()));
        }
        if structures.len() == 1 {
            return Ok(structures[0].clone());
        }
        validate_inputs(structures)?;
        let frame_count = compatible_frame_count(structures)?;
        let mut merger = Merger::new(structures, frame_count);
        merger.append_all()?;
        merger.finish()
    }
}

struct Merger<'a> {
    sources: &'a [Structure],
    frame_count: usize,
    data: StructureData,
    chunks: ChunkBuilder,
    bonds: BondTableBuilder,
    atom_offset: u32,
    residue_offset: u32,
    entity_offset: u32,
    chain_count: u32,
}

impl<'a> Merger<'a> {
    fn new(sources: &'a [Structure], frame_count: usize) -> Self {
        Self {
            sources,
            frame_count,
            data: StructureData::empty(),
            chunks: ChunkBuilder::new(),
            bonds: BondTableBuilder::new(),
            atom_offset: 0,
            residue_offset: 0,
            entity_offset: 0,
            chain_count: 0,
        }
    }

    fn append_all(&mut self) -> Result<(), Vec<Diagnostic>> {
        self.chunks.start_model(0);
        for source in self.sources {
            self.append_entities(source).map_err(single)?;
            self.append_residues(source).map_err(single)?;
            self.append_chains(source).map_err(single)?;
            self.append_atoms(source).map_err(single)?;
            self.append_bonds(source).map_err(single)?;
            self.advance_offsets(source).map_err(single)?;
        }
        Ok(())
    }

    fn append_entities(&mut self, source: &Structure) -> Result<(), Diagnostic> {
        for entity in source.data().topology.entities.iter() {
            let id = remap_symbol(
                source,
                entity_id(source, entity)?,
                &mut self.data.dictionary,
            )?;
            let description = remap_optional(
                source,
                source.data().topology.entities.description(entity),
                &mut self.data.dictionary,
            )?;
            let sequence = source
                .data()
                .topology
                .entities
                .canonical_sequence(entity)
                .iter()
                .map(|symbol| remap_symbol(source, *symbol, &mut self.data.dictionary))
                .collect::<Result<Vec<_>, _>>()?;
            let kind = source
                .data()
                .topology
                .entities
                .kind(entity)
                .ok_or_else(invariant)?;
            self.data
                .topology
                .entities
                .push(id, kind, description, &sequence)
                .map_err(|_| invariant())?;
        }
        Ok(())
    }

    fn append_residues(&mut self, source: &Structure) -> Result<(), Diagnostic> {
        for residue in 0..u32_of(source.residue_count())? {
            let residue = ResidueIndex::new(residue);
            let table = &source.data().topology.residues;
            let atoms = table.atoms(residue).ok_or_else(invariant)?;
            let record = ResidueRecord {
                label_comp_id: remap_symbol(
                    source,
                    table.label_comp_id(residue).ok_or_else(invariant)?,
                    &mut self.data.dictionary,
                )?,
                auth_comp_id: remap_optional(
                    source,
                    table.auth_comp_id(residue),
                    &mut self.data.dictionary,
                )?,
                label_seq_id: table.label_seq_id(residue).into(),
                auth_seq_id: table.auth_seq_id(residue).into(),
                ins_code: remap_optional(
                    source,
                    table.ins_code(residue),
                    &mut self.data.dictionary,
                )?,
                het: table.is_het(residue),
            };
            self.data
                .topology
                .residues
                .push(
                    record,
                    checked_add(atoms.start, self.atom_offset)?
                        ..checked_add(atoms.end, self.atom_offset)?,
                )
                .map_err(|_| invariant())?;
        }
        Ok(())
    }

    fn append_chains(&mut self, source: &Structure) -> Result<(), Diagnostic> {
        for chain in source.data().topology.chains.iter() {
            let table = &source.data().topology.chains;
            let residues = table.residues(chain).ok_or_else(invariant)?;
            let entity = table.entity(chain).ok_or_else(invariant)?;
            let record = ChainRecord {
                label_asym_id: remap_symbol(
                    source,
                    table.label_asym_id(chain).ok_or_else(invariant)?,
                    &mut self.data.dictionary,
                )?,
                auth_asym_id: remap_optional(
                    source,
                    table.auth_asym_id(chain),
                    &mut self.data.dictionary,
                )?,
                entity: EntityIndex::new(checked_add(entity.get(), self.entity_offset)?),
                polymer_kind: table.polymer_kind(chain).ok_or_else(invariant)?,
            };
            self.data
                .topology
                .chains
                .push(
                    record,
                    checked_add(residues.start, self.residue_offset)?
                        ..checked_add(residues.end, self.residue_offset)?,
                )
                .map_err(|_| invariant())?;
            self.chain_count = checked_add(self.chain_count, 1)?;
        }
        Ok(())
    }

    fn append_atoms(&mut self, source: &Structure) -> Result<(), Diagnostic> {
        let positions = source.positions();
        for chunk in source.data().chunks.iter() {
            for old in chunk.atoms() {
                let local = old.checked_sub(chunk.atoms().start).ok_or_else(invariant)?;
                let position = positions.get(old as usize).copied();
                let record = chunk
                    .record(local, &source.data().topology.residues, position)
                    .ok_or_else(invariant)?;
                let remapped = self.remap_atom(source, record)?;
                self.chunks.push(remapped);
            }
        }
        Ok(())
    }

    fn remap_atom(
        &mut self,
        source: &Structure,
        record: AtomRecord,
    ) -> Result<AtomRecord, Diagnostic> {
        Ok(AtomRecord {
            atom_name: remap_symbol(source, record.atom_name, &mut self.data.dictionary)?,
            auth_atom_name: remap_optional(
                source,
                record.auth_atom_name.get(),
                &mut self.data.dictionary,
            )?,
            alternate_component_id: remap_optional(
                source,
                record.alternate_component_id.get(),
                &mut self.data.dictionary,
            )?,
            alt_id: match record.alt_id.symbol() {
                Some(symbol) => {
                    AltId::labelled(remap_symbol(source, symbol, &mut self.data.dictionary)?)
                        .ok_or_else(invariant)?
                }
                None => AltId::BLANK,
            },
            residue: ResidueIndex::new(checked_add(record.residue.get(), self.residue_offset)?),
            ..record
        })
    }

    fn append_bonds(&mut self, source: &Structure) -> Result<(), Diagnostic> {
        for bond in source.data().bonds.iter() {
            self.bonds.push(BondRecord {
                atom_a: AtomIndex::new(checked_add(bond.atom_a.get(), self.atom_offset)?),
                atom_b: AtomIndex::new(checked_add(bond.atom_b.get(), self.atom_offset)?),
                ..bond
            });
        }
        Ok(())
    }

    fn advance_offsets(&mut self, source: &Structure) -> Result<(), Diagnostic> {
        self.atom_offset = checked_add(self.atom_offset, source.atom_count())?;
        self.residue_offset = checked_add(self.residue_offset, u32_of(source.residue_count())?)?;
        self.entity_offset = checked_add(self.entity_offset, u32_of(source.entity_count())?)?;
        Ok(())
    }

    fn finish(mut self) -> Result<Structure, Vec<Diagnostic>> {
        let annotations = merge_annotations(self.sources, &mut self.data.dictionary)?;
        let (chunks, first) = self.chunks.finish();
        self.data.entry = EntryMetadata::default();
        self.data.chunks = chunks.into();
        self.data.bonds = self.bonds.finish_with_availability(
            self.sources
                .iter()
                .all(|source| source.data().bonds.is_available()),
        );
        self.data.annotations = annotations;
        self.data.coords = merged_coordinates(self.sources, self.frame_count, first)?;
        self.data.cell = common_cell(self.sources);
        self.data.generation = CoordinateGeneration::INITIAL;
        add_models(
            &mut self.data.topology,
            self.sources,
            self.frame_count,
            self.chain_count,
        )?;
        let findings = validate(&self.data);
        if findings.is_empty() {
            Ok(Structure::new(self.data))
        } else {
            Err(findings)
        }
    }
}

fn validate_inputs(structures: &[Structure]) -> Result<(), Vec<Diagnostic>> {
    let mut findings: Vec<_> = structures
        .iter()
        .enumerate()
        .flat_map(|(source, structure)| {
            validate(structure.data())
                .into_iter()
                .map(move |finding| finding.with_context("merge source", source.to_string()))
        })
        .collect();
    findings.extend(
        structures
            .iter()
            .enumerate()
            .filter(|(_, structure)| !structure.extensions().is_empty())
            .map(|(source, structure)| {
                Diagnostic::new(Code::E3014)
                    .with_context("merge source", source.to_string())
                    .with_context(
                        "extensions",
                        structure.extensions().keys().collect::<Vec<_>>().join(","),
                    )
            }),
    );
    if findings.is_empty() {
        Ok(())
    } else {
        Err(findings)
    }
}

fn compatible_frame_count(structures: &[Structure]) -> Result<usize, Vec<Diagnostic>> {
    let expected = structures[0].model_count();
    let compatible = structures.iter().all(|structure| {
        !matches!(structure.data().coords, CoordinateStore::Ragged { .. })
            && structure.model_count() == expected
    });
    if compatible {
        Ok(expected)
    } else {
        Err(single(Diagnostic::new(Code::E3012)))
    }
}

fn add_models(
    topology: &mut Topology,
    sources: &[Structure],
    count: usize,
    chains: u32,
) -> Result<(), Vec<Diagnostic>> {
    for model in 0..count {
        let index = ModelIndex::new(u32_of(model).map_err(single)?);
        let stored = sources[0].data().topology.models.model_num(index);
        let Some(number) = stored else {
            return Err(single(invariant()));
        };
        topology
            .models
            .push(number, 0..chains)
            .map_err(|_| single(invariant()))?;
    }
    Ok(())
}

fn merged_coordinates(
    sources: &[Structure],
    frame_count: usize,
    first: CoordinateBlock,
) -> Result<CoordinateStore, Vec<Diagnostic>> {
    if frame_count == 1 {
        return Ok(CoordinateStore::Single(first));
    }
    let mut frames = Vec::with_capacity(frame_count);
    frames.push(first);
    for frame in 1..frame_count {
        let mut merged = CoordinateBlock::new();
        for source in sources {
            let index = ModelIndex::new(u32_of(frame).map_err(single)?);
            if let Some(positions) = source.model_positions(index) {
                for position in positions {
                    merged.push(*position);
                }
            }
        }
        frames.push(merged);
    }
    Ok(CoordinateStore::Dense { frames })
}

fn u32_of(value: usize) -> Result<u32, Diagnostic> {
    u32::try_from(value).map_err(|_| invariant())
}

fn common_cell(sources: &[Structure]) -> Option<UnitCell> {
    let first = sources.first()?.data().cell?;
    sources
        .iter()
        .all(|source| source.data().cell == Some(first))
        .then_some(first)
}

fn entity_id(source: &Structure, entity: EntityIndex) -> Result<SymbolId, Diagnostic> {
    source
        .data()
        .topology
        .entities
        .id(entity)
        .ok_or_else(invariant)
}

fn remap_optional(
    source: &Structure,
    symbol: Option<SymbolId>,
    dictionary: &mut Interner,
) -> Result<OptionalSymbol, Diagnostic> {
    symbol
        .map(|symbol| remap_symbol(source, symbol, dictionary))
        .transpose()
        .map(OptionalSymbol::from)
}

pub(super) fn remap_symbol(
    source: &Structure,
    symbol: SymbolId,
    dictionary: &mut Interner,
) -> Result<SymbolId, Diagnostic> {
    let text = source.resolve(symbol).ok_or_else(invariant)?;
    dictionary.intern(text).map_err(|_| {
        Diagnostic::new(Code::E1901)
            .with_context("limit", "dictionary entries")
            .with_context("identifier", text)
    })
}

fn checked_add(left: u32, right: u32) -> Result<u32, Diagnostic> {
    left.checked_add(right)
        .ok_or_else(|| Diagnostic::new(Code::E1901).with_context("limit", "u32 table rows"))
}

fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}

pub(super) fn single(finding: Diagnostic) -> Vec<Diagnostic> {
    vec![finding]
}

#[cfg(test)]
#[path = "merge_tests.rs"]
mod tests;
