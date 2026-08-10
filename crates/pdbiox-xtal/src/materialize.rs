//! Explicit expansion of a lazy biological assembly.

use crate::AssemblyView;
use pdbiox_core::annotation::{AnnotationColumn, AtomAnnotation, AtomAnnotations};
use pdbiox_core::bond::{BondRecord, BondTableBuilder};
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::structure::CoordinateStore;
use pdbiox_core::topology::{ChainRecord, ResidueRecord};
use pdbiox_core::{
    AtomIndex, ChainIndex, Code, CoordinateBlock, CoordinateGeneration, Diagnostic, ModelIndex,
    ResidueIndex, Structure, StructureData,
};
use pdbiox_geom::Rigid;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Range;

/// Per-atom annotation holding the generated chain's stable instance identifier.
pub const INSTANCE_ID_ANNOTATION: &str = "pdbiox.instance_id";

#[derive(Clone, Debug)]
struct CopySpan {
    source_chain: ChainIndex,
    transform: Rigid,
    source_atoms: Range<u32>,
    output_start: u32,
}

impl AssemblyView {
    /// Explicitly copies topology and transformed coordinates into a structure.
    ///
    /// The source snapshot is untouched. Generated chains receive stable,
    /// collision-free labels and every atom carries [`INSTANCE_ID_ANNOTATION`].
    ///
    /// # Errors
    ///
    /// Returns diagnostics for ragged models, representational overflow, or a
    /// violated source invariant.
    pub fn materialize(&self) -> Result<Structure, Vec<Diagnostic>> {
        if matches!(self.structure.data().coords, CoordinateStore::Ragged { .. }) {
            return Err(single(Diagnostic::new(Code::E6003)));
        }
        Materializer::new(self).build()
    }
}

struct Materializer<'a> {
    view: &'a AssemblyView,
    data: StructureData,
    chunks: ChunkBuilder,
    copies: Vec<CopySpan>,
    source_atoms: Vec<u32>,
    instance_ids: Vec<i64>,
    used_labels: BTreeSet<Box<str>>,
    label_counts: BTreeMap<Box<str>, u32>,
}

impl<'a> Materializer<'a> {
    fn new(view: &'a AssemblyView) -> Self {
        let source = view.source().data();
        let mut data = StructureData::empty();
        data.entry = source.entry.clone();
        data.dictionary = source.dictionary.clone();
        data.cell = source.cell;
        data.topology.entities = source.topology.entities.clone();
        Self {
            view,
            data,
            chunks: ChunkBuilder::new(),
            copies: Vec::with_capacity(view.instance_count()),
            source_atoms: Vec::new(),
            instance_ids: Vec::new(),
            used_labels: BTreeSet::new(),
            label_counts: BTreeMap::new(),
        }
    }

    fn build(mut self) -> Result<Structure, Vec<Diagnostic>> {
        self.chunks.start_model(0);
        for instance in self.view.instances.iter() {
            self.append_instance(instance).map_err(single)?;
        }
        let (chunks, first_frame) = std::mem::take(&mut self.chunks).finish();
        self.data.chunks = chunks.into();
        self.data.coords = self.coordinates(first_frame).map_err(single)?;
        self.data.annotations = self.annotations().map_err(single)?;
        self.data.bonds = self.bonds().map_err(single)?;
        self.data.generation = CoordinateGeneration::INITIAL;
        self.add_models().map_err(single)?;
        let findings = pdbiox_core::structure::validate(&self.data);
        if findings.is_empty() {
            Ok(Structure::new(self.data))
        } else {
            Err(findings)
        }
    }

    fn append_instance(
        &mut self,
        instance: &crate::view::InstanceRecord,
    ) -> Result<(), Diagnostic> {
        let source = self.view.source();
        let source_chains = &source.data().topology.chains;
        let residues = source_chains
            .residues(instance.source_chain)
            .ok_or_else(invariant)?;
        let output_residue_start =
            u32::try_from(self.data.topology.residues.len()).map_err(|_| capacity("residues"))?;
        let output_atom_start =
            u32::try_from(self.source_atoms.len()).map_err(|_| capacity("atoms"))?;
        for source_residue in residues.clone() {
            self.append_residue(ResidueIndex::new(source_residue))?;
        }
        let output_residue_end =
            u32::try_from(self.data.topology.residues.len()).map_err(|_| capacity("residues"))?;
        let label = source_chains
            .label_asym_id(instance.source_chain)
            .and_then(|symbol| source.resolve(symbol))
            .ok_or_else(invariant)?;
        let generated = self.unique_label(label)?;
        let symbol = self
            .data
            .dictionary
            .intern(&generated)
            .map_err(|_| capacity("dictionary entries"))?;
        let entity = source_chains
            .entity(instance.source_chain)
            .ok_or_else(invariant)?;
        let polymer_kind = source_chains
            .polymer_kind(instance.source_chain)
            .ok_or_else(invariant)?;
        self.data
            .topology
            .chains
            .push(
                ChainRecord {
                    label_asym_id: symbol,
                    auth_asym_id: OptionalSymbol::some(symbol),
                    entity,
                    polymer_kind,
                },
                output_residue_start..output_residue_end,
            )
            .map_err(|error| capacity("chains").with_context("cause", error.to_string()))?;
        let output_atom_end =
            u32::try_from(self.source_atoms.len()).map_err(|_| capacity("atoms"))?;
        self.copies.push(CopySpan {
            source_chain: instance.source_chain,
            transform: self.view.transforms[instance.transform],
            source_atoms: instance.atoms.clone(),
            output_start: output_atom_start,
        });
        let instance_atom_count =
            u32::try_from(instance.atoms.len()).map_err(|_| capacity("atoms"))?;
        if output_atom_end
            .checked_sub(output_atom_start)
            .ok_or_else(invariant)?
            != instance_atom_count
        {
            return Err(invariant());
        }
        Ok(())
    }

    fn append_residue(&mut self, source_residue: ResidueIndex) -> Result<(), Diagnostic> {
        let source = self.view.source();
        let table = &source.data().topology.residues;
        let source_atoms = table.atoms(source_residue).ok_or_else(invariant)?;
        let output_residue = ResidueIndex::new(
            u32::try_from(self.data.topology.residues.len()).map_err(|_| capacity("residues"))?,
        );
        let output_start = u32::try_from(self.source_atoms.len()).map_err(|_| capacity("atoms"))?;
        for source_atom in source_atoms.clone() {
            let mut record = atom_record(source, AtomIndex::new(source_atom))?;
            record.position = record
                .position
                .map(|position| self.current_transform().apply(position));
            record.residue = output_residue;
            let local = source_atom
                .checked_sub(source_atoms.start)
                .ok_or_else(invariant)?;
            record.atom_site_id = output_start
                .checked_add(local)
                .and_then(|value| value.checked_add(1))
                .ok_or_else(|| capacity("atoms"))?;
            self.chunks.push(record);
            self.source_atoms.push(source_atom);
            self.instance_ids.push(self.current_instance_id());
        }
        let output_end = u32::try_from(self.source_atoms.len()).map_err(|_| capacity("atoms"))?;
        let record = ResidueRecord {
            label_comp_id: table.label_comp_id(source_residue).ok_or_else(invariant)?,
            auth_comp_id: optional(table.auth_comp_id(source_residue)),
            label_seq_id: table.label_seq_id(source_residue).into(),
            auth_seq_id: table.auth_seq_id(source_residue).into(),
            ins_code: optional(table.ins_code(source_residue)),
            het: table.is_het(source_residue),
        };
        self.data
            .topology
            .residues
            .push(record, output_start..output_end)
            .map_err(|error| capacity("residues").with_context("cause", error.to_string()))?;
        Ok(())
    }

    fn current_transform(&self) -> Rigid {
        let position = self.copies.len();
        let instance = &self.view.instances[position];
        self.view.transforms[instance.transform]
    }

    fn current_instance_id(&self) -> i64 {
        i64::from(self.view.instances[self.copies.len()].instance_id.get())
    }

    fn unique_label(&mut self, source: &str) -> Result<Box<str>, Diagnostic> {
        let count = self.label_counts.entry(source.into()).or_insert(0);
        loop {
            *count = count
                .checked_add(1)
                .ok_or_else(|| capacity("chain identifiers"))?;
            let candidate: Box<str> = if *count == 1 {
                source.into()
            } else {
                format!("{source}-{count}").into_boxed_str()
            };
            if self.used_labels.insert(candidate.clone()) {
                return Ok(candidate);
            }
        }
    }

    fn coordinates(&self, first: CoordinateBlock) -> Result<CoordinateStore, Diagnostic> {
        let count = self.view.source().model_count();
        if count == 1 {
            return Ok(CoordinateStore::Single(first));
        }
        let mut frames = Vec::with_capacity(count);
        frames.push(first);
        for model in 1..count {
            let source = self
                .view
                .source()
                .model_positions(ModelIndex::new(
                    u32::try_from(model).map_err(|_| capacity("models"))?,
                ))
                .ok_or_else(invariant)?;
            let mut frame = CoordinateBlock::with_capacity(self.source_atoms.len());
            for copy in &self.copies {
                for atom in copy.source_atoms.clone() {
                    let position = source.get(atom as usize).copied().ok_or_else(invariant)?;
                    frame.push(copy.transform.apply(position));
                }
            }
            frames.push(frame);
        }
        Ok(CoordinateStore::Dense { frames })
    }

    fn annotations(&self) -> Result<AtomAnnotations, Diagnostic> {
        let mut output = AtomAnnotations::default();
        for (name, annotation) in self.view.source().annotations().iter() {
            let copied = copy_annotation(annotation, &self.source_atoms)?;
            let _ = output.insert(name, copied);
        }
        let instance_ids = AnnotationColumn::from_values(self.instance_ids.clone())
            .map_err(|error| invariant().with_context("cause", error.to_string()))?;
        let _ = output.insert(
            INSTANCE_ID_ANNOTATION,
            AtomAnnotation::Integer(instance_ids),
        );
        Ok(output)
    }

    fn bonds(&self) -> Result<pdbiox_core::BondTable, Diagnostic> {
        let source = self.view.source();
        let mut builder = BondTableBuilder::new();
        let groups = copy_groups(&self.copies);
        for bond in source.data().bonds.iter() {
            let chain_a = atom_chain(source, bond.atom_a)?;
            let chain_b = atom_chain(source, bond.atom_b)?;
            for group in groups.values() {
                let Some(left) = group.get(&chain_a.get()) else {
                    continue;
                };
                let Some(right) = group.get(&chain_b.get()) else {
                    continue;
                };
                if chain_a == chain_b {
                    for copy in left {
                        builder.push(remap_bond(bond, &self.copies[*copy], &self.copies[*copy])?);
                    }
                } else {
                    for copy_a in left {
                        for copy_b in right {
                            builder.push(remap_bond(
                                bond,
                                &self.copies[*copy_a],
                                &self.copies[*copy_b],
                            )?);
                        }
                    }
                }
            }
        }
        Ok(builder.finish_with_availability(source.data().bonds.is_available()))
    }

    fn add_models(&mut self) -> Result<(), Diagnostic> {
        let chains =
            u32::try_from(self.data.topology.chains.len()).map_err(|_| capacity("chains"))?;
        for model in 0..self.view.source().model_count() {
            let index = ModelIndex::new(u32::try_from(model).map_err(|_| capacity("models"))?);
            let number = self
                .view
                .source()
                .data()
                .topology
                .models
                .model_num(index)
                .ok_or_else(|| missing_model_number(model))?;
            self.data
                .topology
                .models
                .push(number, 0..chains)
                .map_err(|error| capacity("models").with_context("cause", error.to_string()))?;
        }
        Ok(())
    }
}

fn missing_model_number(model: usize) -> Diagnostic {
    Diagnostic::new(Code::E4105)
        .with_context("model", model.to_string())
        .with_context("required", "explicit source model number")
}

fn atom_record(source: &Structure, atom: AtomIndex) -> Result<AtomRecord, Diagnostic> {
    let candidate = source
        .data()
        .chunks
        .partition_point(|chunk| chunk.atoms().end <= atom.get());
    let chunk = source.data().chunks.get(candidate).ok_or_else(invariant)?;
    let local = atom
        .get()
        .checked_sub(chunk.atoms().start)
        .ok_or_else(invariant)?;
    let position = source.positions().get(atom.as_usize()).copied();
    chunk
        .record(local, &source.data().topology.residues, position)
        .ok_or_else(invariant)
}

fn copy_annotation(
    annotation: &AtomAnnotation,
    source_atoms: &[u32],
) -> Result<AtomAnnotation, Diagnostic> {
    match annotation {
        AtomAnnotation::Boolean(column) => {
            copy_column(column, source_atoms).map(AtomAnnotation::Boolean)
        }
        AtomAnnotation::Integer(column) => {
            copy_column(column, source_atoms).map(AtomAnnotation::Integer)
        }
        AtomAnnotation::Real(column) => copy_column(column, source_atoms).map(AtomAnnotation::Real),
        AtomAnnotation::Symbol(column) => {
            copy_column(column, source_atoms).map(AtomAnnotation::Symbol)
        }
        _ => {
            Err(Diagnostic::new(Code::E3013)
                .with_context("annotation", "unsupported physical type"))
        }
    }
}

fn copy_column<T: Copy>(
    column: &AnnotationColumn<T>,
    source_atoms: &[u32],
) -> Result<AnnotationColumn<T>, Diagnostic> {
    let entries = source_atoms
        .iter()
        .map(|atom| column.get(*atom).ok_or_else(invariant))
        .collect::<Result<Vec<_>, _>>()?;
    AnnotationColumn::from_entries(entries)
        .map_err(|error| invariant().with_context("cause", error.to_string()))
}

type TransformKey = [u64; 12];

fn copy_groups(copies: &[CopySpan]) -> BTreeMap<TransformKey, BTreeMap<u32, Vec<usize>>> {
    let mut groups = BTreeMap::<TransformKey, BTreeMap<u32, Vec<usize>>>::new();
    for (position, copy) in copies.iter().enumerate() {
        groups
            .entry(transform_key(copy.transform))
            .or_default()
            .entry(copy.source_chain.get())
            .or_default()
            .push(position);
    }
    groups
}

fn transform_key(transform: Rigid) -> TransformKey {
    let mut key = [0_u64; 12];
    for row in 0..3 {
        for column in 0..3 {
            key[row * 3 + column] = transform.rotation[row][column].to_bits();
        }
        key[9 + row] = transform.translation[row].to_bits();
    }
    key
}

fn atom_chain(source: &Structure, atom: AtomIndex) -> Result<ChainIndex, Diagnostic> {
    let residue = source
        .data()
        .topology
        .residues
        .containing(atom.get())
        .ok_or_else(invariant)?;
    source
        .data()
        .topology
        .chains
        .containing(residue.get())
        .ok_or_else(invariant)
}

fn remap_bond(
    bond: BondRecord,
    left: &CopySpan,
    right: &CopySpan,
) -> Result<BondRecord, Diagnostic> {
    Ok(BondRecord {
        atom_a: remap_atom(bond.atom_a, left)?,
        atom_b: remap_atom(bond.atom_b, right)?,
        ..bond
    })
}

fn remap_atom(atom: AtomIndex, copy: &CopySpan) -> Result<AtomIndex, Diagnostic> {
    if !copy.source_atoms.contains(&atom.get()) {
        return Err(invariant());
    }
    let local = atom
        .get()
        .checked_sub(copy.source_atoms.start)
        .ok_or_else(invariant)?;
    copy.output_start
        .checked_add(local)
        .map(AtomIndex::new)
        .ok_or_else(|| capacity("atoms"))
}

fn optional(symbol: Option<pdbiox_core::SymbolId>) -> OptionalSymbol {
    match symbol {
        Some(symbol) => OptionalSymbol::some(symbol),
        None => OptionalSymbol::NONE,
    }
}

fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}

fn capacity(name: &str) -> Diagnostic {
    Diagnostic::new(Code::E1901).with_context("limit", name)
}

fn single(finding: Diagnostic) -> Vec<Diagnostic> {
    vec![finding]
}

#[cfg(test)]
#[path = "materialize_tests.rs"]
mod tests;
