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
        if matches!(
            &self.structure.data().coords,
            CoordinateStore::Ragged { .. }
        ) {
            return Err(single(Diagnostic::new(Code::E6003)));
        }

        let materializer = Materializer::new(self).map_err(single)?;
        materializer.build()
    }
}

/// Transactional builder for an explicitly materialized assembly.
struct Materializer<'a> {
    view: &'a AssemblyView,
    data: StructureData,
    chunks: ChunkBuilder,
    copies: Vec<CopySpan>,
    source_atoms: Vec<u32>,
    instance_ids: Vec<i64>,
    used_labels: IdentityHashSet<SymbolId>,
    label_counts: IdentityHashMap<SymbolId, u32>,
}

struct ResidueAppend<'a> {
    source: &'a Structure,
    transform: Rigid,
    instance_id: i64,
    expected_source_atom: &'a mut u32,
    output_atom: &'a mut u32,
    chunk_cursor: &'a mut usize,
}

impl<'a> Materializer<'a> {
    /// Creates a materializer and reserves storage for the expected expansion.
    ///
    /// Source metadata that remains valid after expansion is shared or cloned
    /// into the independent output snapshot.
    fn new(view: &'a AssemblyView) -> Result<Self, Diagnostic> {
        let source = view.source().data();
        let total_atoms = expanded_atom_count(view)?;

        let mut data = StructureData::empty();
        data.entry = source.entry.clone();
        data.dictionary = source.dictionary.clone();
        data.cell = source.cell;
        data.topology.entities = source.topology.entities.clone();

        let mut copies = Vec::new();
        copies
            .try_reserve_exact(view.instance_count())
            .map_err(|_| capacity("assembly instances"))?;

        let mut source_atoms = Vec::new();
        source_atoms
            .try_reserve_exact(total_atoms)
            .map_err(|_| capacity("atoms"))?;

        let mut instance_ids = Vec::new();
        instance_ids
            .try_reserve_exact(total_atoms)
            .map_err(|_| capacity("atoms"))?;

        let mut used_labels = IdentityHashSet::default();
        used_labels
            .try_reserve(view.instance_count())
            .map_err(|_| capacity("chain identifiers"))?;

        let mut label_counts = IdentityHashMap::default();
        label_counts
            .try_reserve(view.instance_count())
            .map_err(|_| capacity("chain identifiers"))?;

        Ok(Self {
            view,
            data,
            chunks: ChunkBuilder::new(),
            copies,
            source_atoms,
            instance_ids,
            used_labels,
            label_counts,
        })
    }

    /// Builds and validates the fully independent assembly structure.
    fn build(mut self) -> Result<Structure, Vec<Diagnostic>> {
        self.chunks.start_model(0);

        for instance in self.view.instances.as_ref() {
            self.append_instance(instance).map_err(single)?;
        }

        let (chunks, first_frame) = std::mem::take(&mut self.chunks).finish();
        self.data.chunks = chunks.into();

        let coordinates = self.coordinates(first_frame).map_err(single)?;
        let annotations = self.annotations().map_err(single)?;
        let bonds = self.bonds().map_err(single)?;

        self.data.coords = coordinates;
        self.data.annotations = annotations;
        self.data.bonds = bonds;
        self.data.generation = CoordinateGeneration::INITIAL;

        self.add_models().map_err(single)?;

        let findings = pdbiox_core::structure::validate(&self.data);

        if findings.is_empty() {
            Ok(Structure::new(self.data))
        } else {
            Err(findings)
        }
    }

    /// Appends one transformed source-chain instance to the output hierarchy.
    ///
    /// The declared instance atom span is checked against the source hierarchy
    /// while atoms are copied, preventing equal-length but misaligned ranges
    /// from producing chemically inconsistent later coordinate frames.
    fn append_instance(
        &mut self,
        instance: &crate::view::InstanceRecord,
    ) -> Result<(), Diagnostic> {
        if instance.atoms.start > instance.atoms.end {
            return Err(invariant());
        }

        let view = self.view;
        let source = view.source();
        let source_chains = &source.data().topology.chains;

        let transform = view
            .transforms
            .get(instance.transform)
            .copied()
            .ok_or_else(invariant)?;
        let instance_id = i64::from(instance.instance_id.get());

        let residues = source_chains
            .residues(instance.source_chain)
            .ok_or_else(invariant)?;

        let output_residue_start =
            u32::try_from(self.data.topology.residues.len()).map_err(|_| capacity("residues"))?;
        let output_atom_start =
            u32::try_from(self.source_atoms.len()).map_err(|_| capacity("atoms"))?;

        let mut output_atom = output_atom_start;
        let mut expected_source_atom = instance.atoms.start;
        let mut chunk_cursor = source
            .data()
            .chunks
            .partition_point(|chunk| chunk.atoms().end <= instance.atoms.start);
        let mut append = ResidueAppend {
            source,
            transform,
            instance_id,
            expected_source_atom: &mut expected_source_atom,
            output_atom: &mut output_atom,
            chunk_cursor: &mut chunk_cursor,
        };

        for source_residue in residues {
            self.append_residue(ResidueIndex::new(source_residue), &mut append)?;
        }

        if expected_source_atom != instance.atoms.end {
            return Err(invariant());
        }

        let output_atom_end =
            u32::try_from(self.source_atoms.len()).map_err(|_| capacity("atoms"))?;

        if output_atom != output_atom_end {
            return Err(invariant());
        }

        let output_residue_end =
            u32::try_from(self.data.topology.residues.len()).map_err(|_| capacity("residues"))?;

        let source_symbol = source_chains
            .label_asym_id(instance.source_chain)
            .ok_or_else(invariant)?;
        let source_label = source.resolve(source_symbol).ok_or_else(invariant)?;
        let generated_symbol = self.unique_label(source_symbol, source_label)?;

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
                    label_asym_id: generated_symbol,
                    auth_asym_id: OptionalSymbol::some(generated_symbol),
                    entity,
                    polymer_kind,
                },
                output_residue_start..output_residue_end,
            )
            .map_err(|error| capacity("chains").with_context("cause", error.to_string()))?;

        self.copies.push(CopySpan {
            source_chain: instance.source_chain,
            transform,
            source_atoms: instance.atoms.clone(),
            output_start: output_atom_start,
        });

        Ok(())
    }

    /// Copies one source residue and all of its atoms into the current instance.
    ///
    /// Transform, instance identifier, output index, and chunk cursor are
    /// supplied by the caller so they are not rediscovered for every atom.
    fn append_residue(
        &mut self,
        source_residue: ResidueIndex,
        append: &mut ResidueAppend<'_>,
    ) -> Result<(), Diagnostic> {
        let table = &append.source.data().topology.residues;
        let source_range = table.atoms(source_residue).ok_or_else(invariant)?;

        let output_residue = ResidueIndex::new(
            u32::try_from(self.data.topology.residues.len()).map_err(|_| capacity("residues"))?,
        );

        let output_start = *append.output_atom;

        for source_atom in source_range {
            if source_atom != *append.expected_source_atom {
                return Err(invariant());
            }

            let atom = AtomIndex::new(source_atom);
            let mut record = atom_record(append.source, atom, append.chunk_cursor)?;

            record.position = record
                .position
                .map(|position| append.transform.apply(position));
            record.residue = output_residue;

            let atom_site_id = append
                .output_atom
                .checked_add(1)
                .ok_or_else(|| capacity("atoms"))?;
            record.atom_site_id = atom_site_id;

            self.chunks.push(record);
            self.source_atoms.push(source_atom);
            self.instance_ids.push(append.instance_id);

            *append.output_atom = atom_site_id;
            *append.expected_source_atom = append
                .expected_source_atom
                .checked_add(1)
                .ok_or_else(|| capacity("atoms"))?;
        }

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
            .push(record, output_start..*append.output_atom)
            .map_err(|error| capacity("residues").with_context("cause", error.to_string()))?;

        Ok(())
    }

    /// Generates a stable collision-free label and returns its dictionary symbol.
    ///
    /// Existing source symbols are reused for the first available label.
    /// Generated suffixes are interned only when collisions require them.
    fn unique_label(
        &mut self,
        source_symbol: SymbolId,
        source: &str,
    ) -> Result<SymbolId, Diagnostic> {
        let mut count = match self.label_counts.get(&source_symbol) {
            Some(count) => *count,
            None => 0,
        };

        loop {
            count = count
                .checked_add(1)
                .ok_or_else(|| capacity("chain identifiers"))?;

            let symbol = if count == 1 {
                source_symbol
            } else {
                let candidate = format!("{source}-{count}");
                self.data
                    .dictionary
                    .intern(&candidate)
                    .map_err(|_| capacity("dictionary entries"))?
            };

            if self.used_labels.insert(symbol) {
                self.label_counts.insert(source_symbol, count);
                return Ok(symbol);
            }
        }
    }

    /// Materializes transformed coordinates for every source model.
    ///
    /// The first frame is already produced by the chunk builder. Remaining
    /// frames replay the recorded copy spans in output atom order.
    fn coordinates(&self, first: CoordinateBlock) -> Result<CoordinateStore, Diagnostic> {
        let source = self.view.source();
        let model_count = source.model_count();

        if model_count == 0 {
            return Err(invariant());
        }

        if model_count == 1 {
            return Ok(CoordinateStore::Single(first));
        }

        let mut frames = Vec::new();
        frames
            .try_reserve_exact(model_count)
            .map_err(|_| capacity("models"))?;
        frames.push(first);

        for model in 1..model_count {
            let model = u32::try_from(model).map_err(|_| capacity("models"))?;
            let source_positions = source
                .model_positions(ModelIndex::new(model))
                .ok_or_else(invariant)?;

            let mut frame = CoordinateBlock::with_capacity(self.source_atoms.len());

            for copy in &self.copies {
                for atom in copy.source_atoms.clone() {
                    let index = usize::try_from(atom).map_err(|_| capacity("atoms"))?;
                    let position = source_positions.get(index).copied().ok_or_else(invariant)?;

                    frame.push(copy.transform.apply(position));
                }
            }

            frames.push(frame);
        }

        Ok(CoordinateStore::Dense { frames })
    }

    /// Copies source annotation columns and installs generated instance IDs.
    ///
    /// The instance-ID buffer is moved into its annotation column rather than
    /// cloned after assembly expansion.
    fn annotations(&mut self) -> Result<AtomAnnotations, Diagnostic> {
        let mut output = AtomAnnotations::default();

        for (name, annotation) in self.view.source().annotations().iter() {
            let copied = copy_annotation(annotation, &self.source_atoms)?;
            let _ = output.insert(name, copied);
        }

        let instance_ids = std::mem::take(&mut self.instance_ids);
        let instance_ids = AnnotationColumn::from_values(instance_ids)
            .map_err(|error| invariant().with_context("cause", error.to_string()))?;

        let _ = output.insert(
            INSTANCE_ID_ANNOTATION,
            AtomAnnotation::Integer(instance_ids),
        );

        Ok(output)
    }

    /// Reconstructs source bonds between assembly copies sharing a transform.
    ///
    /// Transform-group matches are cached per source-chain pair so repeated
    /// bonds between the same chains do not rescan every assembly transform.
    fn bonds(&self) -> Result<pdbiox_core::BondTable, Diagnostic> {
        let source = self.view.source();
        let groups = copy_groups(&self.copies);
        let mut matching_groups = BTreeMap::<(u32, u32), Vec<usize>>::new();
        let mut builder = BondTableBuilder::new();

        for bond in source.data().bonds.iter() {
            let chain_a = atom_chain(source, bond.atom_a)?;
            let chain_b = atom_chain(source, bond.atom_b)?;
            let key = chain_pair_key(chain_a, chain_b);

            let group_indices = matching_groups
                .entry(key)
                .or_insert_with(|| groups_for_pair(&groups, chain_a.get(), chain_b.get()));

            for &group_index in group_indices.iter() {
                let group = groups.get(group_index).ok_or_else(invariant)?;
                let left = group.get(&chain_a.get()).ok_or_else(invariant)?;
                let right = group.get(&chain_b.get()).ok_or_else(invariant)?;

                if chain_a == chain_b {
                    for &copy in left {
                        let span = self.copies.get(copy).ok_or_else(invariant)?;
                        builder.push(remap_bond(bond, span, span)?);
                    }
                } else {
                    for &copy_a in left {
                        let left_span = self.copies.get(copy_a).ok_or_else(invariant)?;

                        for &copy_b in right {
                            let right_span = self.copies.get(copy_b).ok_or_else(invariant)?;
                            builder.push(remap_bond(bond, left_span, right_span)?);
                        }
                    }
                }
            }
        }

        Ok(builder.finish_with_availability(source.data().bonds.is_available()))
    }

    /// Copies explicit source model numbers onto the materialized topology.
    fn add_models(&mut self) -> Result<(), Diagnostic> {
        let view = self.view;
        let source = view.source();

        let chains =
            u32::try_from(self.data.topology.chains.len()).map_err(|_| capacity("chains"))?;

        for model in 0..source.model_count() {
            let index = ModelIndex::new(u32::try_from(model).map_err(|_| capacity("models"))?);

            let number = source
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

