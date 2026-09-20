//! Large, deterministic Arrow materialisation and streaming resource probes.

use super::{ResourceRecord, measure_case};
use arrow::ffi_stream::ArrowArrayStreamReader;

const LARGE_ARROW_ATOMS: u32 = 2_097_152;
const LARGE_ARROW_RESIDUE_ATOMS: u32 = 16;

pub(super) fn run_arrow_stream_large() -> Result<ResourceRecord, String> {
    let structure = large_arrow_structure()?;
    let table = molframe::interop::AtomTable::new(structure.engine());
    measure_case("arrow_stream_large", || {
        let stream = table
            .arrow_stream()
            .map_err(|error| format!("Arrow stream creation failed: {error}"))?;
        let mut reader = ArrowArrayStreamReader::try_new(stream.into_ffi())
            .map_err(|error| format!("Arrow stream import failed: {error}"))?;
        reader.try_fold(0_u64, |rows, batch| {
            let batch = batch.map_err(|error| format!("Arrow stream pull failed: {error}"))?;
            add_rows(rows, batch.num_rows())
        })
    })
}

pub(super) fn run_arrow_batches_large() -> Result<ResourceRecord, String> {
    let structure = large_arrow_structure()?;
    let table = molframe::interop::AtomTable::new(structure.engine());
    measure_case("arrow_batches_large", || {
        let batches = table
            .record_batches()
            .map_err(|error| format!("Arrow materialisation failed: {error}"))?;
        batches
            .iter()
            .try_fold(0_u64, |rows, batch| add_rows(rows, batch.num_rows()))
    })
}

fn add_rows(rows: u64, batch_rows: usize) -> Result<u64, String> {
    let batch_rows =
        u64::try_from(batch_rows).map_err(|_| "Arrow batch row count exceeds u64".to_owned())?;
    rows.checked_add(batch_rows)
        .ok_or_else(|| "Arrow row count exceeds u64".to_owned())
}

fn large_arrow_structure() -> Result<molframe::Structure, String> {
    use molframe_core::optional::{OptionalI32, OptionalSymbol};
    use molframe_core::topology::ResidueRecord;

    let mut data = molframe::engine::core::StructureData::empty();
    let atom_name = data
        .dictionary
        .intern("CA")
        .map_err(|error| format!("atom name interning failed: {error}"))?;
    let component = data
        .dictionary
        .intern("ALA")
        .map_err(|error| format!("component interning failed: {error}"))?;
    let mut chunks = molframe::engine::core::ChunkBuilder::new();
    let atom_capacity = usize::try_from(LARGE_ARROW_ATOMS)
        .map_err(|_| "Arrow stress atom count exceeds usize".to_owned())?;
    chunks.reserve(atom_capacity);
    let mut first = 0_u32;
    while first < LARGE_ARROW_ATOMS {
        let end = first
            .checked_add(LARGE_ARROW_RESIDUE_ATOMS)
            .map_or(LARGE_ARROW_ATOMS, |end| end.min(LARGE_ARROW_ATOMS));
        let sequence = i32::try_from(first / LARGE_ARROW_RESIDUE_ATOMS)
            .map_err(|_| "Arrow stress residue index exceeds i32".to_owned())?
            .checked_add(1)
            .ok_or_else(|| "Arrow stress residue index exceeds i32".to_owned())?;
        let residue = data
            .topology
            .residues
            .push(
                ResidueRecord {
                    label_comp_id: component,
                    auth_comp_id: OptionalSymbol::NONE,
                    label_seq_id: OptionalI32::some(sequence),
                    auth_seq_id: OptionalI32::some(sequence),
                    ins_code: OptionalSymbol::NONE,
                    het: false,
                },
                first..end,
            )
            .map_err(|error| format!("Arrow stress residue failed: {error}"))?;
        for atom in first..end {
            chunks.push(molframe::engine::core::AtomRecord {
                position: Some([0.0, 0.0, 0.0]),
                element: molframe::Element::CARBON,
                atom_name,
                auth_atom_name: OptionalSymbol::NONE,
                alternate_component_id: OptionalSymbol::NONE,
                alt_id: molframe::AltId::BLANK,
                residue,
                occupancy: (1.0, molframe::engine::core::Presence::Present),
                b_factor: (10.0, molframe::engine::core::Presence::Present),
                formal_charge: (0, molframe::engine::core::Presence::Inapplicable),
                atom_site_id: atom + 1,
            });
        }
        first = end;
    }
    let (atom_chunks, coordinates) = chunks.finish();
    data.chunks = atom_chunks.into();
    data.coords = molframe::engine::core::CoordinateStore::Single(coordinates);
    Ok(molframe::Structure::from(
        molframe_core::structure::Structure::new(data),
    ))
}
