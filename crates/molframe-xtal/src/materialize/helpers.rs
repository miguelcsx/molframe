/// Computes the expected expanded atom count with overflow checking.
fn expanded_atom_count(view: &AssemblyView) -> Result<usize, Diagnostic> {
    view.instances.iter().try_fold(0usize, |total, instance| {
        total
            .checked_add(instance.atoms.len())
            .ok_or_else(|| capacity("atoms"))
    })
}

/// Creates a diagnostic for a missing explicit source model number.
fn missing_model_number(model: usize) -> Diagnostic {
    Diagnostic::new(Code::E4105)
        .with_context("model", model.to_string())
        .with_context("required", "explicit source model number")
}

/// Reconstructs one source atom record using a monotonic chunk cursor.
///
/// Assembly instance atoms are copied in source order, allowing chunk lookup to
/// advance linearly instead of binary-searching the chunk table for every atom.
fn atom_record(
    source: &Structure,
    atom: AtomIndex,
    chunk_cursor: &mut usize,
) -> Result<AtomRecord, Diagnostic> {
    let chunks = &source.data().chunks;

    loop {
        let chunk = chunks.get(*chunk_cursor).ok_or_else(invariant)?;
        let atoms = chunk.atoms();

        if atom.get() < atoms.start {
            return Err(invariant());
        }

        if atom.get() < atoms.end {
            let local = atom.get().checked_sub(atoms.start).ok_or_else(invariant)?;
            let position = source.positions().get(atom.as_usize()).copied();

            return chunk
                .record(local, &source.data().topology.residues, position)
                .ok_or_else(invariant);
        }

        *chunk_cursor = chunk_cursor
            .checked_add(1)
            .ok_or_else(|| capacity("chunks"))?;
    }
}

/// Copies one annotation column according to the output-to-source atom mapping.
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

/// Copies entries from one physical annotation column in output atom order.
fn copy_column<T: Copy>(
    column: &AnnotationColumn<T>,
    source_atoms: &[u32],
) -> Result<AnnotationColumn<T>, Diagnostic> {
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(source_atoms.len())
        .map_err(|_| capacity("annotation rows"))?;

    for &atom in source_atoms {
        entries.push(column.get(atom).ok_or_else(invariant)?);
    }

    AnnotationColumn::from_entries(entries)
        .map_err(|error| invariant().with_context("cause", error.to_string()))
}

/// Groups assembly copies sharing an exactly identical rigid transform.
fn copy_groups(copies: &[CopySpan]) -> Vec<CopyGroup> {
    let mut grouped = BTreeMap::<TransformKey, CopyGroup>::new();

    for (position, copy) in copies.iter().enumerate() {
        grouped
            .entry(transform_key(copy.transform))
            .or_default()
            .entry(copy.source_chain.get())
            .or_default()
            .push(position);
    }

    grouped.into_values().collect()
}

/// Finds transform groups containing both source chains of a bond.
fn groups_for_pair(groups: &[CopyGroup], chain_a: u32, chain_b: u32) -> Vec<usize> {
    groups
        .iter()
        .enumerate()
        .filter_map(|(index, group)| {
            (group.contains_key(&chain_a) && group.contains_key(&chain_b)).then_some(index)
        })
        .collect()
}

/// Produces an orientation-independent cache key for a source-chain pair.
fn chain_pair_key(left: ChainIndex, right: ChainIndex) -> (u32, u32) {
    let left = left.get();
    let right = right.get();

    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

/// Produces an exact deterministic key for a rigid transform.
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

/// Resolves the source chain owning one atom.
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

/// Remaps both endpoints of one source bond into materialized copies.
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

/// Remaps one source atom into the output span of a copied chain.
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

/// Converts an optional symbol into its compact topology representation.
fn optional(symbol: Option<SymbolId>) -> OptionalSymbol {
    match symbol {
        Some(symbol) => OptionalSymbol::some(symbol),
        None => OptionalSymbol::NONE,
    }
}

/// Creates the internal-invariant diagnostic used for impossible assembly state.
fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}

/// Creates a representational-capacity diagnostic for the named resource.
fn capacity(name: &str) -> Diagnostic {
    Diagnostic::new(Code::E1901).with_context("limit", name)
}

/// Wraps one diagnostic for APIs that report a collection of findings.
fn single(finding: Diagnostic) -> Vec<Diagnostic> {
    vec![finding]
}

