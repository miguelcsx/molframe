//! Criterion coverage for element tables and CCD lowering.

use criterion::{Criterion, black_box};
use molframe_bench::{ccd_atp, ccd_hem, input};
use molframe_chem::{
    ComponentProvider, PeoeAtom, PeoeAtomType, PeoeBond, PeoeOptions, RadiusSet, SmartsPattern,
    component_peoe_charges, element_properties, equivalence_classes, peoe_charges, read_ccd,
    vdw_radius,
};
use molframe_core::contract::DictionaryVersion;

fn bench_element_table(c: &mut Criterion) {
    c.bench_function("chem_element_radius_table", |b| {
        b.iter(|| {
            for element in [
                molframe_core::Element::CARBON,
                molframe_core::Element::NITROGEN,
                molframe_core::Element::OXYGEN,
                molframe_core::Element::SULFUR,
            ] {
                black_box((
                    element_properties(element),
                    vdw_radius(element, RadiusSet::Bondi),
                ));
            }
        });
    });
}

fn bench_ccd(c: &mut Criterion) {
    let buffer = input(ccd_hem());
    c.bench_function("chem_read_ccd/HEM", |b| {
        b.iter(|| black_box(read_ccd(&buffer, DictionaryVersion::new("bench"))));
    });
}

fn bench_component_kernels(c: &mut Criterion) {
    let buffer = input(ccd_atp());
    let provider = match read_ccd(&buffer, DictionaryVersion::new("bench")) {
        Ok((provider, _)) => provider,
        Err(findings) => panic!("ATP benchmark dictionary failed: {findings:?}"),
    };
    let component = match provider.get("ATP") {
        Ok(Some(component)) => component,
        Ok(None) => panic!("ATP benchmark dictionary omitted ATP"),
        Err(error) => panic!("ATP benchmark lookup failed: {error}"),
    };
    if let Err(error) = component_peoe_charges(&component, PeoeOptions::default()) {
        panic!("ATP PEOE benchmark setup failed: {error}");
    }
    let carbon_path = match SmartsPattern::parse("C-C-O") {
        Ok(pattern) => pattern,
        Err(error) => panic!("SMARTS benchmark pattern failed: {error}"),
    };
    let mut group = c.benchmark_group("chem_component_kernels");
    group.bench_function("equivalence_classes/ATP", |b| {
        b.iter(|| black_box(equivalence_classes(&component)));
    });
    group.bench_function("peoe/ATP", |b| {
        b.iter(|| black_box(component_peoe_charges(&component, PeoeOptions::default())));
    });
    group.bench_function("smarts_exists/ATP", |b| {
        b.iter(|| black_box(carbon_path.matches(&component)));
    });
    group.bench_function("smarts_find_all/ATP", |b| {
        b.iter(|| black_box(carbon_path.find_matches(&component)));
    });
    group.finish();

    let atoms = vec![
        PeoeAtom {
            atom_type: PeoeAtomType::CSp3,
            formal_charge: 0.0,
        };
        4_096
    ];
    let bonds: Vec<PeoeBond> = (1..atoms.len())
        .map(|atom| PeoeBond {
            atom_a: atom - 1,
            atom_b: atom,
        })
        .collect();
    c.bench_function("chem_peoe/synthetic_chain_4096", |b| {
        b.iter(|| black_box(peoe_charges(&atoms, &bonds, PeoeOptions::default())));
    });
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_element_table(&mut criterion);
    bench_ccd(&mut criterion);
    bench_component_kernels(&mut criterion);
    criterion.final_summary();
}
