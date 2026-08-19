//! Process-isolated allocation and memory measurements for representative workflows.

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicU64, Ordering};

use serde::Serialize;

struct CountingAllocator;

static ALLOCATION_COUNT: AtomicU64 = AtomicU64::new(0);
static ALLOCATED_BYTES: AtomicU64 = AtomicU64::new(0);
static CURRENT_BYTES: AtomicU64 = AtomicU64::new(0);
static PEAK_BYTES: AtomicU64 = AtomicU64::new(0);

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

// SAFETY: every operation delegates to the platform allocator. The counters
// are independent atomic bookkeeping and do not change allocation semantics.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        // SAFETY: the caller supplies the layout required by GlobalAlloc.
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            record_allocation(layout.size());
        }
        pointer
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        // SAFETY: the caller supplies the layout required by GlobalAlloc.
        let pointer = unsafe { System.alloc_zeroed(layout) };
        if !pointer.is_null() {
            record_allocation(layout.size());
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        // SAFETY: the pointer and layout are the pair previously returned by
        // this allocator, as required by GlobalAlloc.
        unsafe { System.dealloc(pointer, layout) };
        subtract_current(layout.size());
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        // SAFETY: the caller supplies the live allocation and its old layout.
        let replacement = unsafe { System.realloc(pointer, layout, size) };
        if !replacement.is_null() {
            ALLOCATION_COUNT.fetch_add(1, Ordering::SeqCst);
            ALLOCATED_BYTES.fetch_add(size as u64, Ordering::SeqCst);
            if size >= layout.size() {
                add_current(size - layout.size());
            } else {
                subtract_current(layout.size() - size);
            }
        }
        replacement
    }
}

fn record_allocation(size: usize) {
    ALLOCATION_COUNT.fetch_add(1, Ordering::SeqCst);
    ALLOCATED_BYTES.fetch_add(size as u64, Ordering::SeqCst);
    add_current(size);
}

fn add_current(size: usize) {
    let current = CURRENT_BYTES.fetch_add(size as u64, Ordering::SeqCst) + size as u64;
    let mut peak = PEAK_BYTES.load(Ordering::SeqCst);
    while current > peak {
        match PEAK_BYTES.compare_exchange(peak, current, Ordering::SeqCst, Ordering::SeqCst) {
            Ok(_) => break,
            Err(observed) => peak = observed,
        }
    }
}

fn subtract_current(size: usize) {
    CURRENT_BYTES
        .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |current| {
            Some(current.saturating_sub(size as u64))
        })
        .ok();
}

fn signed_difference(after: u64, before: u64) -> i64 {
    if after >= before {
        match i64::try_from(after - before) {
            Ok(value) => value,
            Err(_) => i64::MAX,
        }
    } else {
        match i64::try_from(before - after) {
            Ok(value) => -value,
            Err(_) => i64::MIN,
        }
    }
}

#[derive(Clone, Copy)]
struct CounterSnapshot {
    allocations: u64,
    allocated_bytes: u64,
    current_bytes: u64,
    peak_bytes: u64,
}

fn snapshot() -> CounterSnapshot {
    CounterSnapshot {
        allocations: ALLOCATION_COUNT.load(Ordering::SeqCst),
        allocated_bytes: ALLOCATED_BYTES.load(Ordering::SeqCst),
        current_bytes: CURRENT_BYTES.load(Ordering::SeqCst),
        peak_bytes: PEAK_BYTES.load(Ordering::SeqCst),
    }
}

#[derive(Serialize)]
struct ResourceRecord {
    schema_version: u32,
    case: &'static str,
    allocation_count: u64,
    allocated_bytes: u64,
    peak_live_bytes: u64,
    live_bytes_delta: i64,
    result_digest: u64,
}

fn measure_case<F>(name: &'static str, operation: F) -> Result<ResourceRecord, String>
where
    F: FnOnce() -> Result<u64, String>,
{
    let before = snapshot();
    PEAK_BYTES.store(before.current_bytes, Ordering::SeqCst);
    let result = black_box(operation());
    let after = snapshot();
    let digest = result?;
    Ok(ResourceRecord {
        schema_version: 1,
        case: name,
        allocation_count: after.allocations.saturating_sub(before.allocations),
        allocated_bytes: after.allocated_bytes.saturating_sub(before.allocated_bytes),
        peak_live_bytes: after.peak_bytes.saturating_sub(before.current_bytes),
        live_bytes_delta: signed_difference(after.current_bytes, before.current_bytes),
        result_digest: digest,
    })
}

fn read_mmcif(bytes: &[u8], name: &'static str) -> Result<pdbiox::Structure, String> {
    pdbiox::read_bytes(bytes.to_vec(), Some(name), &pdbiox::ReadOptions::new())
        .map(|(structure, _)| structure)
        .map_err(|findings| format!("{name} read failed: {findings:?}"))
}

fn run(name: &str) -> Result<ResourceRecord, String> {
    match name {
        "mmcif_read_medium" => measure_case("mmcif_read_medium", || {
            let bytes = pdbiox_bench::Sample::Medium
                .cif()
                .ok_or_else(|| "medium mmCIF fixture is missing".to_owned())?;
            let structure = read_mmcif(bytes, "medium.cif")?;
            Ok(black_box(u64::from(structure.atom_count())))
        }),
        "bcif_read_medium" => measure_case("bcif_read_medium", || {
            let bytes = pdbiox_bench::Sample::Medium.bcif();
            let structure = read_mmcif(bytes, "medium.bcif")?;
            Ok(black_box(u64::from(structure.atom_count())))
        }),
        "cif_lower_medium" => measure_case("cif_lower_medium", || {
            let bytes = pdbiox_bench::Sample::Medium
                .cif()
                .ok_or_else(|| "medium mmCIF fixture is missing".to_owned())?;
            let input = pdbiox::InputBuffer::from_bytes(bytes.to_vec());
            let document = pdbiox::cif::parse(&input)
                .map_err(|findings| format!("CIF parse failed: {findings:?}"))?
                .0;
            let structure = pdbiox::cif::lower(&document, &pdbiox::ReadOptions::new())
                .map_err(|findings| format!("CIF lowering failed: {findings:?}"))?
                .0;
            Ok(black_box(u64::from(structure.atom_count())))
        }),
        "selection_medium" => {
            use pdbiox::QueryStructure as _;

            let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Medium);
            measure_case("selection_medium", || {
                let evaluation = structure
                    .select_text("within 5 of element H", &pdbiox::AnalysisPolicy::default())
                    .map_err(|findings| format!("selection failed: {findings:?}"))?;
                Ok(black_box(evaluation.selection.len() as u64))
            })
        }
        "contacts_medium" => {
            let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Medium);
            measure_case("contacts_medium", || {
                let contacts =
                    pdbiox::analysis::atom_contacts(&structure, 4.0, pdbiox::SpatialBackend::Auto)
                        .map_err(|error| format!("contacts failed: {error}"))?;
                Ok(black_box(contacts.len() as u64))
            })
        }
        "sasa_medium" => {
            let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Medium);
            let positions = structure.positions().to_vec();
            let radii = vec![1.7_f32; positions.len()];
            measure_case("sasa_medium", || {
                let areas = pdbiox::surface::shrake_rupley(&positions, &radii, 1.4, 96)
                    .map_err(|error| format!("SASA failed: {error}"))?;
                Ok(black_box(areas.len() as u64))
            })
        }
        "rmsd_medium" => {
            let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Medium);
            let reference = structure.positions().to_vec();
            let mobile = pdbiox_bench::perturbed(&reference, 0.001);
            measure_case("rmsd_medium", || {
                let rmsd = pdbiox::rmsd(&mobile, &reference)
                    .map_err(|error| format!("RMSD failed: {error:?}"))?;
                Ok(black_box(rmsd.to_bits()))
            })
        }
        "coordinate_handoff" => {
            let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Medium);
            measure_case("coordinate_handoff", || {
                let positions = structure.positions();
                black_box(positions.as_ptr());
                Ok(black_box(positions.len() as u64))
            })
        }
        "dlpack_coordinates" => {
            let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Medium);
            measure_case("dlpack_coordinates", || {
                let tensor = pdbiox::DlpackTensor::coordinates(&structure)
                    .map_err(|error| format!("DLPack failed: {error}"))?;
                Ok(black_box(u64::from(tensor.cost() as u8)))
            })
        }
        "trajectory_contacts" => run_trajectory_contacts(),
        _ => Err(format!("unknown resource case: {name}")),
    }
}

fn run_trajectory_contacts() -> Result<ResourceRecord, String> {
    let structure = pdbiox_bench::structure(pdbiox_bench::Sample::Tiny);
    let trajectory = pdbiox::traj::Trajectory::from_frames(
        (0_u16..64)
            .map(|frame| pdbiox::traj::Frame {
                positions: structure
                    .positions()
                    .iter()
                    .map(|&[x, y, z]| [x + f32::from(frame) * 1.0e-4, y, z])
                    .collect(),
            })
            .collect(),
    );
    let policy = pdbiox::AnalysisPolicy::default();
    let kernel = pdbiox::analysis::contacts_kernel(3.0, pdbiox::SpatialBackend::Auto);
    measure_case("trajectory_contacts", || {
        let analysis =
            pdbiox::analysis::analyse_trajectory(&structure, &trajectory, &policy, &kernel, 4)
                .map_err(|error| format!("trajectory analysis failed: {error}"))?;
        black_box(analysis.value);
        Ok(64)
    })
}

fn main() -> Result<(), String> {
    let case = std::env::args()
        .nth(1)
        .ok_or_else(|| "usage: pdbiox-resource-bench CASE".to_owned())?;
    let record = run(&case)?;
    println!(
        "{}",
        serde_json::to_string(&record).map_err(|error| format!("resource JSON failed: {error}"))?
    );
    Ok(())
}
