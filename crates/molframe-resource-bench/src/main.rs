//! Process-isolated allocation and memory measurements for representative workflows.

mod arrow_bench;
mod generated_bcif_bench;
mod generated_structure_bench;
mod kernel_bench;
mod modelcif_bench;
mod mrc_bench;
mod stream_bench;
mod xtc_bench;

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::path::Path;
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

fn stabilize_shared_runtime() -> Result<(), String> {
    let context = molframe::ExecutionContext::default();
    let plan = molframe_core::parallel::BlockPlan::new(context.worker_budget(), 1);
    let warmed = molframe_core::parallel::map_blocks_in(plan, &context, |index, _| index)
        .map_err(|error| format!("shared runtime warm-up failed: {error}"))?;
    black_box(warmed);
    black_box(context);
    Ok(())
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
    #[serde(skip_serializing_if = "Option::is_none")]
    model_count: Option<u64>,
}

fn measure_case<F>(name: &'static str, operation: F) -> Result<ResourceRecord, String>
where
    F: FnOnce() -> Result<u64, String>,
{
    stabilize_shared_runtime()?;
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
        model_count: None,
    })
}

fn measure_retained_case<F, T>(name: &'static str, operation: F) -> Result<ResourceRecord, String>
where
    F: FnOnce() -> Result<(u64, Option<u64>, T), String>,
{
    stabilize_shared_runtime()?;
    let before = snapshot();
    PEAK_BYTES.store(before.current_bytes, Ordering::SeqCst);
    let (digest, model_count, retained) = black_box(operation())?;
    black_box(&retained);
    let after = snapshot();
    let record = ResourceRecord {
        schema_version: 1,
        case: name,
        allocation_count: after.allocations.saturating_sub(before.allocations),
        allocated_bytes: after.allocated_bytes.saturating_sub(before.allocated_bytes),
        peak_live_bytes: after.peak_bytes.saturating_sub(before.current_bytes),
        live_bytes_delta: signed_difference(after.current_bytes, before.current_bytes),
        result_digest: digest,
        model_count,
    };
    drop(retained);
    Ok(record)
}

fn read_mmcif(bytes: &[u8], name: &'static str) -> Result<molframe::Structure, String> {
    molframe::read_bytes(bytes.to_vec(), Some(name), &molframe::ReadOptions::new())
        .map(|(structure, _)| structure)
        .map_err(|findings| format!("{name} read failed: {findings:?}"))
}

fn run(name: &str) -> Result<ResourceRecord, String> {
    match name {
        "mmcif_read_medium" => measure_case("mmcif_read_medium", || {
            let bytes = molframe_bench::Sample::Medium
                .cif()
                .ok_or_else(|| "medium mmCIF fixture is missing".to_owned())?;
            let structure = read_mmcif(bytes, "medium.cif")?;
            Ok(black_box(u64::from(structure.atom_count())))
        }),
        "bcif_read_medium" => measure_case("bcif_read_medium", || {
            let bytes = molframe_bench::Sample::Medium.bcif();
            let structure = read_mmcif(bytes, "medium.bcif")?;
            Ok(black_box(u64::from(structure.atom_count())))
        }),
        "cif_lower_medium" => measure_case("cif_lower_medium", || {
            let bytes = molframe_bench::Sample::Medium
                .cif()
                .ok_or_else(|| "medium mmCIF fixture is missing".to_owned())?;
            let input = molframe::InputBuffer::from_bytes(bytes.to_vec());
            let document = molframe::formats::cif::parse(&input)
                .map_err(|findings| format!("CIF parse failed: {findings:?}"))?
                .0;
            let structure = molframe::formats::cif::lower(&document, &molframe::ReadOptions::new())
                .map_err(|findings| format!("CIF lowering failed: {findings:?}"))?
                .0;
            Ok(black_box(u64::from(structure.atom_count())))
        }),
        "selection_medium" => {
            use molframe::QueryStructure as _;

            let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
            measure_case("selection_medium", || {
                let evaluation = structure
                    .select_text(
                        "within 5 of element H",
                        &molframe::AnalysisPolicy::default(),
                    )
                    .map_err(|findings| format!("selection failed: {findings:?}"))?;
                Ok(black_box(evaluation.selection.len() as u64))
            })
        }
        "contacts_medium" => {
            let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
            measure_case("contacts_medium", || {
                let contacts = molframe::analysis::atom_contacts(
                    &structure,
                    4.0,
                    molframe::spatial::SpatialBackend::Auto,
                    &molframe::ExecutionContext::default(),
                )
                .map_err(|error| format!("contacts failed: {error}"))?;
                Ok(black_box(contacts.len() as u64))
            })
        }
        "sasa_medium" => {
            let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
            let positions = structure.positions().to_vec();
            let radii = vec![1.7_f32; positions.len()];
            measure_case("sasa_medium", || {
                let areas = molframe::surface::shrake_rupley(
                    &positions,
                    &radii,
                    1.4,
                    96,
                    &molframe::ExecutionContext::default(),
                )
                .map_err(|error| format!("SASA failed: {error}"))?;
                Ok(black_box(areas.len() as u64))
            })
        }
        "rmsd_medium" => {
            let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
            let reference = structure.positions().to_vec();
            let mobile = molframe_bench::perturbed(&reference, 0.001);
            measure_case("rmsd_medium", || {
                let rmsd = molframe::geometry::rmsd(&mobile, &reference)
                    .map_err(|error| format!("RMSD failed: {error:?}"))?;
                Ok(black_box(rmsd.to_bits()))
            })
        }
        "coordinate_handoff" => {
            let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
            measure_case("coordinate_handoff", || {
                let positions = structure.positions();
                black_box(positions.as_ptr());
                Ok(black_box(positions.len() as u64))
            })
        }
        "dlpack_coordinates" => {
            let structure = molframe_bench::structure(molframe_bench::Sample::Medium);
            measure_case("dlpack_coordinates", || {
                let tensor = molframe::interop::DlpackTensor::coordinates(&structure)
                    .map_err(|error| format!("DLPack failed: {error}"))?;
                Ok(black_box(u64::from(tensor.cost() as u8)))
            })
        }
        _ => run_extended(name),
    }
}

fn run_extended(name: &str) -> Result<ResourceRecord, String> {
    match name {
        "trajectory_contacts" => kernel_bench::run_trajectory_contacts(),
        "pore_profile_100000" => kernel_bench::run_pore_profile_100000(),
        "msa_eight_512" => kernel_bench::run_msa_eight_512(),
        "msa_sixty_four_128" => kernel_bench::run_msa_sixty_four_128(),
        "arrow_stream_large" => arrow_bench::run_arrow_stream_large(),
        "arrow_batches_large" => arrow_bench::run_arrow_batches_large(),
        "mrc_block_1g" => mrc_bench::run_mrc_block_1g(),
        "stream_synthetic_256m" => stream_bench::run_stream_synthetic_256m(),
        "stream_synthetic_1g" => stream_bench::run_stream_synthetic_1g(),
        "lex_synthetic_1g" => stream_bench::run_lex_synthetic_1g(),
        "scan_synthetic_256m" => stream_bench::run_scan_synthetic_256m(),
        "scan_synthetic_1g" => stream_bench::run_scan_synthetic_1g(),
        "file_copied_1g" => stream_bench::run_file_copied_1g(),
        "file_mapped_1g" => stream_bench::run_file_mapped_1g(),
        "bcif_stream_file" => measure_case("bcif_stream_file", || {
            let path = std::env::var_os("MOLFRAME_BENCH_FILE")
                .map(std::path::PathBuf::from)
                .ok_or_else(|| "MOLFRAME_BENCH_FILE is not set".to_owned())?;
            stream_bench::drain_structure_batches("bcif_stream_file", &path)
        }),
        _ => Err(format!("unknown resource case: {name}")),
    }
}

fn read_file(
    name: &'static str,
    path: &Path,
    options: &molframe::ReadOptions,
) -> Result<ResourceRecord, String> {
    measure_retained_case(name, || {
        let structure = molframe::read_with_options(path, options)
            .map(|(structure, _)| structure)
            .map_err(|findings| format!("{} read failed: {findings:?}", path.display()))?;
        let model_count = u64::try_from(structure.model_count())
            .map_err(|_| "model count exceeds u64".to_owned())?;
        let atom_rows = total_atom_rows(&structure)?;
        Ok((atom_rows, Some(model_count), structure))
    })
}

fn total_atom_rows(structure: &molframe::Structure) -> Result<u64, String> {
    if let Some(models) = structure.engine().ragged_models() {
        return models.iter().try_fold(0_u64, |total, model| {
            total
                .checked_add(u64::from(model.atom_count()))
                .ok_or_else(|| "total atom-row count exceeds u64".to_owned())
        });
    }
    let models =
        u64::try_from(structure.model_count()).map_err(|_| "model count exceeds u64".to_owned())?;
    u64::from(structure.atom_count())
        .checked_mul(models)
        .ok_or_else(|| "total atom-row count exceeds u64".to_owned())
}

fn require_no_more_arguments(
    arguments: &mut impl Iterator<Item = String>,
    case: &str,
) -> Result<(), String> {
    if arguments.next().is_some() {
        return Err(format!("{case} accepts exactly one path"));
    }
    Ok(())
}

fn write_bcif_file(path: &Path) -> Result<ResourceRecord, String> {
    let structure = molframe::read_with_options(path, &molframe::ReadOptions::new())
        .map(|(structure, _)| structure)
        .map_err(|findings| format!("{} read failed: {findings:?}", path.display()))?;
    let atom_count = u64::from(structure.atom_count());
    let model_count =
        u64::try_from(structure.model_count()).map_err(|_| "model count exceeds u64".to_owned())?;
    let options = molframe::formats::cif::CifWriteOptions::new()
        .with_generated_connection_ids()
        .with_connection_type_id("covale");
    measure_retained_case("file_write_bcif", || {
        let bytes = molframe::write_bcif_with_options(&structure, &options)
            .map_err(|findings| format!("BinaryCIF write failed: {findings:?}"))?;
        black_box(bytes.len());
        Ok((atom_count, Some(model_count), bytes))
    })
}

fn main() -> Result<(), String> {
    let mut arguments = std::env::args().skip(1);
    let case = arguments
        .next()
        .ok_or_else(|| "usage: molframe-resource-bench CASE".to_owned())?;
    let record = if case == "generated_bcif_batches" {
        let bytes = arguments
            .next()
            .ok_or_else(|| format!("usage: molframe-resource-bench {case} LOGICAL_BYTES"))?
            .parse::<u64>()
            .map_err(|error| format!("{case}: invalid logical byte count: {error}"))?;
        if arguments.next().is_some() {
            return Err(format!("{case} accepts exactly one logical byte count"));
        }
        generated_bcif_bench::run(bytes)?
    } else if let Some(format) = generated_structure_bench::GeneratedFormat::parse(&case) {
        let bytes = arguments
            .next()
            .ok_or_else(|| format!("usage: molframe-resource-bench {case} BYTES"))?
            .parse::<u64>()
            .map_err(|error| format!("{case}: invalid byte count: {error}"))?;
        if arguments.next().is_some() {
            return Err(format!("{case} accepts exactly one byte count"));
        }
        generated_structure_bench::run(format, bytes)?
    } else if matches!(
        case.as_str(),
        "file_read"
            | "file_read_atomic"
            | "file_write_bcif"
            | "modelcif_file"
            | "mrc_block_file"
            | "structure_batch_file"
            | "window_scan_file"
    ) || case == "xtc_file"
    {
        let path = arguments
            .next()
            .ok_or_else(|| format!("usage: molframe-resource-bench {case} PATH"))?;
        if case == "structure_batch_file" {
            let maximum_spill_bytes = arguments
                .next()
                .ok_or_else(|| {
                    "usage: molframe-resource-bench structure_batch_file PATH MAX_SPILL_BYTES"
                        .to_owned()
                })?
                .parse::<u64>()
                .map_err(|error| format!("invalid spill byte ceiling: {error}"))?;
            if arguments.next().is_some() {
                return Err("structure_batch_file accepts a path and spill ceiling".to_owned());
            }
            stream_bench::run_structure_batch_file(Path::new(&path), maximum_spill_bytes)?
        } else if case == "window_scan_file" {
            require_no_more_arguments(&mut arguments, &case)?;
            stream_bench::run_window_scan_file(Path::new(&path))?
        } else if case == "xtc_file" {
            require_no_more_arguments(&mut arguments, &case)?;
            xtc_bench::read_file(Path::new(&path))?
        } else if case == "modelcif_file" {
            require_no_more_arguments(&mut arguments, &case)?;
            modelcif_bench::read_file(Path::new(&path))?
        } else if case == "mrc_block_file" {
            require_no_more_arguments(&mut arguments, &case)?;
            mrc_bench::read_file("mrc_block_file", Path::new(&path))?
        } else if case == "file_write_bcif" {
            require_no_more_arguments(&mut arguments, &case)?;
            write_bcif_file(Path::new(&path))?
        } else {
            require_no_more_arguments(&mut arguments, &case)?;
            let options =
                molframe::ReadOptions::new().only_atomic_coords(case == "file_read_atomic");
            read_file(
                if case == "file_read" {
                    "file_read"
                } else {
                    "file_read_atomic"
                },
                Path::new(&path),
                &options,
            )?
        }
    } else {
        if arguments.next().is_some() {
            return Err(format!(
                "resource case {case} does not accept extra arguments"
            ));
        }
        run(&case)?
    };
    println!(
        "{}",
        serde_json::to_string(&record).map_err(|error| format!("resource JSON failed: {error}"))?
    );
    Ok(())
}
