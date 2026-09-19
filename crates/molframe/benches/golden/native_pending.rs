use criterion::BenchmarkGroup;

#[path = "native_pending/analysis.rs"]
mod analysis;
#[path = "native_pending/compare_interop.rs"]
mod compare_interop;
#[path = "native_pending/trajectory.rs"]
mod trajectory;
#[path = "native_pending/xtal.rs"]
mod xtal;

pub(super) fn register(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    analysis::register(group);
    compare_interop::register(group);
    trajectory::register(group);
    xtal::register(group);
}
