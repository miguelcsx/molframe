//! Criterion coverage for bounded policy-space expansion.

use std::collections::BTreeSet;

use criterion::{Criterion, black_box};
use pdbiox_audit::{PolicyDimension, PolicySpace, audit_batch};
use pdbiox_core::contract::{AnalysisPolicy, MissingPolicy, Namespace};

fn bench_policy_plan(c: &mut Criterion) {
    let space = PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::identifiers([
            Namespace::Auth,
            Namespace::Label,
        ]))
        .vary(PolicyDimension::identifiers([
            Namespace::Auth,
            Namespace::Label,
        ]));
    c.bench_function("audit_policy_plan", |b| {
        b.iter(|| black_box(space.clone().plan()));
    });

    let plan = match PolicySpace::new(AnalysisPolicy::default())
        .vary(PolicyDimension::missing_atoms([
            MissingPolicy::Report,
            MissingPolicy::Indeterminate,
        ]))
        .plan()
    {
        Ok(plan) => plan,
        Err(error) => panic!("audit benchmark plan failed: {error}"),
    };
    let subjects: Vec<u32> = (0..1_024).collect();
    c.bench_function("audit_batch/1024_subjects", |b| {
        b.iter(|| {
            black_box(audit_batch(
                &plan,
                &subjects,
                |subject: &u32, policy: &AnalysisPolicy| {
                    let mut values = BTreeSet::from([*subject]);
                    if policy.missing_atoms == MissingPolicy::Report {
                        values.insert(subject.saturating_add(1));
                    }
                    Ok::<_, ()>(values)
                },
                Clone::clone,
            ))
        });
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_policy_plan(&mut criterion);
    criterion.final_summary();
}
