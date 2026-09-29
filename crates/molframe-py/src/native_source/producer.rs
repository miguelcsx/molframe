//! Builds a structure's capsule and answers its C callbacks.

use super::abi::{
    ABI_VERSION, CAPSULE_NAME, NativeAtom, NativeBond, NativeCapsule, NativeSourceV2,
    NativeTopology,
};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use std::ffi::c_char;
use std::slice;
use std::sync::{Mutex, OnceLock};

pub(crate) fn capsule<'py>(
    py: Python<'py>,
    structure: &molframe::Structure,
) -> PyResult<Bound<'py, PyCapsule>> {
    let topology = compact_topology(structure);
    let coordinates = structure.coordinates();
    let api = NativeSourceV2 {
        abi_version: ABI_VERSION,
        coordinate_generation: structure.engine().generation().get(),
        coordinates: coordinates.as_ptr().cast::<f32>(),
        coordinate_rows: coordinates.len(),
        atoms: topology.atoms.as_ptr(),
        atom_count: topology.atoms.len(),
        residue_atom_start: topology.residue_atom_start.as_ptr(),
        residue_atom_start_len: topology.residue_atom_start.len(),
        chain_residue_start: topology.chain_residue_start.as_ptr(),
        chain_residue_start_len: topology.chain_residue_start.len(),
        model_chain_start: topology.model_chain_start.as_ptr(),
        model_chain_start_len: topology.model_chain_start.len(),
        bonds: topology.bonds.as_ptr(),
        bond_count: topology.bonds.len(),
        select,
        encode_bcif,
    };
    PyCapsule::new_with_value(
        py,
        NativeCapsule {
            api,
            structure: structure.clone(),
            topology,
            encoded_bcif: OnceLock::new(),
            queries: Mutex::new(crate::query_cache::QueryCache::default()),
            pending: Mutex::new(crate::query_cache::PendingRows::default()),
        },
        CAPSULE_NAME,
    )
}

/// The canonical-write decisions the capsule payload makes on its own.
///
/// The payload is a transport between two extensions, not a file anyone keeps,
/// so it may invent what canonical output otherwise refuses to: a block
/// identifier for a structure read from a format with no `_entry.id`, and
/// `_struct_conn` identifiers and a generic connection type for the bond graph,
/// which retains neither. Without them a structure that carries bonds, or no
/// entry identifier, could not cross the boundary at all.
fn transport_options() -> molframe::formats::cif::CifWriteOptions {
    molframe::formats::cif::CifWriteOptions::new()
        .with_block_id("molframe")
        .with_generated_connection_ids()
        .with_connection_type_id("covale")
}

unsafe extern "C" fn encode_bcif(
    source: *const NativeSourceV2,
    output: *mut u8,
    capacity: usize,
) -> i64 {
    // SAFETY: the API is the first field of its retained capsule allocation.
    let capsule = unsafe { &*source.cast::<NativeCapsule>() };
    let encoded = capsule.encoded_bcif.get_or_init(|| {
        molframe::write_bcif_with_options(&capsule.structure, &transport_options()).map_err(|_| ())
    });
    let Ok(encoded) = encoded else {
        return -1;
    };
    let Ok(length) = i64::try_from(encoded.len()) else {
        return -1;
    };
    if output.is_null() || capacity == 0 {
        return length;
    }
    if capacity < encoded.len() {
        return -1;
    }
    // SAFETY: the caller supplies at least `encoded.len()` writable bytes.
    unsafe { std::ptr::copy_nonoverlapping(encoded.as_ptr(), output, encoded.len()) };
    length
}

fn compact_topology(structure: &molframe::Structure) -> NativeTopology {
    let data = structure.engine().data();
    let atoms = data
        .atoms()
        .map(|atom| NativeAtom {
            element: atom
                .element()
                .map_or(0, |element| u16::from(element.atomic_number())),
            residue: atom.residue().map_or(0, |residue| residue.index().get()),
        })
        .collect();
    let residue_atom_start = offsets(
        data.residues()
            .map(|residue| residue.atoms().map(|atom| atom.index().get())),
    );
    let chain_residue_start = offsets(
        data.chains()
            .map(|chain| chain.residues().map(|residue| residue.index().get())),
    );
    let model_chain_start = offsets(
        data.models()
            .map(|model| model.chains().map(|chain| chain.index().get())),
    );
    let bonds = structure
        .bonds()
        .iter()
        .map(|bond| NativeBond {
            first: bond.atom_a.get(),
            second: bond.atom_b.get(),
            aromatic: u8::from(bond.order == molframe::BondOrder::Aromatic),
            order: bond_order_code(bond.order),
            reserved: [0; 2],
        })
        .collect();
    NativeTopology {
        atoms,
        residue_atom_start,
        chain_residue_start,
        model_chain_start,
        bonds,
    }
}

const fn bond_order_code(order: molframe::BondOrder) -> u8 {
    match order {
        molframe::BondOrder::Unknown => NativeBond::ORDER_UNKNOWN,
        molframe::BondOrder::Single => NativeBond::ORDER_SINGLE,
        molframe::BondOrder::Double => NativeBond::ORDER_DOUBLE,
        molframe::BondOrder::Triple => NativeBond::ORDER_TRIPLE,
        molframe::BondOrder::Quadruple => NativeBond::ORDER_QUADRUPLE,
        molframe::BondOrder::Aromatic => NativeBond::ORDER_AROMATIC,
        molframe::BondOrder::Polymeric => NativeBond::ORDER_POLYMERIC,
    }
}

fn offsets<I, R>(rows: I) -> Vec<u32>
where
    I: Iterator<Item = R>,
    R: Iterator<Item = u32>,
{
    let mut starts = Vec::new();
    let mut end = 0;
    for row in rows {
        let mut row = row.peekable();
        let mut start = end;
        if let Some(value) = row.peek().copied() {
            start = value;
        }
        starts.push(start);
        end = row.last().map_or(start, |value| value.saturating_add(1));
    }
    starts.push(end);
    starts
}

/// The callback's answer to a null argument.
///
/// Producers that predate the size request return this for a null output
/// buffer, which is how a consumer tells them apart.
pub(super) const SELECT_NULL_ARGUMENT: i64 = -1;

/// Evaluates one query for a consumer.
///
/// A null `output` with zero `output_len` is a size request: the query is
/// evaluated, its rows are kept for the copy request that follows, and the
/// row count is returned. Any other call copies rows into `output`, taking the
/// kept rows when they answer the same query and evaluating otherwise. Each
/// distinct query text is compiled once per capsule.
unsafe extern "C" fn select(
    api: *const NativeSourceV2,
    source: *const c_char,
    source_len: usize,
    output: *mut u32,
    output_len: usize,
) -> i64 {
    if api.is_null() || source.is_null() {
        return SELECT_NULL_ARGUMENT;
    }
    let size_request = output.is_null() && output_len == 0;
    if output.is_null() && !size_request {
        return SELECT_NULL_ARGUMENT;
    }
    // SAFETY: the consumer passes back the API pointer supplied by this
    // producer; `NativeSourceV2` is the first field of `NativeCapsule`.
    let capsule = unsafe { &*api.cast::<NativeCapsule>() };
    // SAFETY: the caller promises a readable query buffer of `source_len`.
    let bytes = unsafe { slice::from_raw_parts(source.cast::<u8>(), source_len) };
    let Ok(source) = std::str::from_utf8(bytes) else {
        return -2;
    };
    let kept = match capsule.pending.lock() {
        Ok(mut pending) => pending.take(source),
        Err(_) => None,
    };
    let evaluation = match kept {
        Some(evaluation) => evaluation,
        None => match evaluate(capsule, source) {
            Some(evaluation) => evaluation,
            None => return -3,
        },
    };
    let Ok(count) = usize::try_from(evaluation.selection.len()) else {
        return -4;
    };
    let Ok(encoded) = i64::try_from(count) else {
        return -5;
    };
    if size_request {
        if let Ok(mut pending) = capsule.pending.lock() {
            pending.keep(source, evaluation);
        }
        return encoded;
    }
    if count > output_len {
        return -4;
    }
    // SAFETY: the consumer supplies a writable buffer of `output_len` rows.
    let output = unsafe { slice::from_raw_parts_mut(output, output_len) };
    for (slot, atom) in output.iter_mut().zip(evaluation.selection.iter()) {
        *slot = atom;
    }
    encoded
}

/// Evaluates `source` with a compiled query from the capsule's cache.
fn evaluate(capsule: &NativeCapsule, source: &str) -> Option<molframe::query::Evaluation> {
    use molframe::QueryStructure;

    let policy = molframe::AnalysisPolicy::default();
    let query = match capsule.queries.lock() {
        Ok(mut cache) => cache.get_or_compile(source).ok()?,
        // A poisoned cache only loses reuse; the query is still answered.
        Err(_) => std::sync::Arc::new(molframe::Query::compile(source).ok()?),
    };
    capsule.structure.select_query(&query, &policy).ok()
}
