use super::command::metric_name;
use crate::MetricChoice;

#[test]
fn metric_names_are_stable_machine_vocabulary() {
    assert_eq!(metric_name(MetricChoice::Lddt), "lddt");
    assert_eq!(metric_name(MetricChoice::TmScore), "tm_score");
    assert_eq!(metric_name(MetricChoice::GdtTs), "gdt_ts");
    assert_eq!(metric_name(MetricChoice::GdtHa), "gdt_ha");
    assert_eq!(metric_name(MetricChoice::DockQ), "dockq");
}
