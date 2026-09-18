/// One global MMTF bond translated to model-local atom indices.
#[derive(Clone, Copy)]
struct GlobalBond {
    atom_a: u32,
    atom_b: u32,
    order: BondOrder,
}

/// Partitions global bonds once instead of scanning the complete bond table for
/// every model. Cross-model bonds retain the previous behavior and are ignored.
fn partition_global_bonds(
    decoded: &mut Decoded,
    ranges: &[ModelRange],
) -> Result<Vec<Vec<GlobalBond>>, Diagnostic> {
    if !decoded.global_bonds.len().is_multiple_of(2) {
        return Err(schema_error("global bond endpoint list has odd length"));
    }
    let mut counts = vec![0usize; ranges.len()];
    for pair in decoded.global_bonds.as_chunks::<2>().0 {
        if let Some((model, _, _)) = local_bond_endpoints(pair, ranges)? {
            counts[model] = counts[model]
                .checked_add(1)
                .ok_or_else(|| schema_error("MMTF model bond count overflows"))?;
        }
    }
    let mut partitioned = counts
        .into_iter()
        .map(Vec::with_capacity)
        .collect::<Vec<_>>();
    for (bond_index, pair) in decoded.global_bonds.as_chunks::<2>().0.iter().enumerate() {
        let Some((model, atom_a, atom_b)) = local_bond_endpoints(pair, ranges)? else {
            continue;
        };
        partitioned[model].push(GlobalBond {
            atom_a,
            atom_b,
            order: decoded_order(
                &decoded.global_orders,
                &decoded.global_resonance,
                bond_index,
            ),
        });
    }
    decoded.global_bonds = Vec::new();
    decoded.global_orders = Vec::new();
    decoded.global_resonance = Vec::new();
    Ok(partitioned)
}

/// Resolves both global endpoints and returns their common model and offsets.
fn local_bond_endpoints(
    pair: &[i32],
    ranges: &[ModelRange],
) -> Result<Option<(usize, u32, u32)>, Diagnostic> {
    let first = usize_of(pair[0])?;
    let second = usize_of(pair[1])?;
    let Some((first_model, first_local)) = model_atom(first, ranges)? else {
        return Ok(None);
    };
    let Some((second_model, second_local)) = model_atom(second, ranges)? else {
        return Ok(None);
    };
    if first_model == second_model {
        Ok(Some((first_model, first_local, second_local)))
    } else {
        Ok(None)
    }
}

/// Binary-searches the disjoint, ascending model atom ranges.
fn model_atom(
    atom: usize,
    ranges: &[ModelRange],
) -> Result<Option<(usize, u32)>, Diagnostic> {
    let Some(model) = ranges
        .partition_point(|range| range.atom_start <= atom)
        .checked_sub(1)
    else {
        return Ok(None);
    };
    let range = ranges[model];
    let end = range
        .atom_start
        .checked_add(range.atom_count)
        .ok_or_else(|| schema_error("MMTF model atom range overflows"))?;
    if atom >= end {
        return Ok(None);
    }
    let local = atom - range.atom_start;
    let local = u32::try_from(local)
        .map_err(|_| schema_error("MMTF local atom index exceeds u32"))?;
    Ok(Some((model, local)))
}
