//! Criterion coverage for superposition-free and fitted structure scores.

use criterion::{Criterion, black_box};
use pdbiox_bench::{Sample, coordinates, perturbed, structure};
use pdbiox_chem::{Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider};
use pdbiox_compare::{
    CeOptions, ComparisonAlignment, ContactArea, DockQOptions, PointMapping, PointMatch, QsOptions,
    cad_score, ce_align, contact_map_similarity, dockq, equivalent_atom_mappings, gdt_ha, gdt_ts,
    interface_rmsd, lddt, ligand_symmetry_rmsd, map_sequence_to_structure, measure_mapping,
    pocket_rmsd, qs_score, tm_score, weighted_rmsd,
};
use pdbiox_core::contract::{DictionaryVersion, Namespace};
use pdbiox_core::{BondOrder, Element};
use std::fmt::Debug;
use std::sync::Arc;

trait BenchRequired<T> {
    fn required(self, context: &str) -> T;
}

impl<T, E: Debug> BenchRequired<T> for Result<T, E> {
    fn required(self, context: &str) -> T {
        match self {
            Ok(value) => value,
            Err(error) => panic!("{context}: {error:?}"),
        }
    }
}

fn bench_scores(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let reference = coordinates(&structure);
    let model = perturbed(&reference, 0.02);
    let mut group = c.benchmark_group("compare_scores");
    group.bench_function("lddt", |b| {
        b.iter(|| black_box(lddt(&model, &reference, 15.0)));
    });
    group.bench_function("tm_score", |b| {
        b.iter(|| black_box(tm_score(&model, &reference)));
    });
    group.bench_function("gdt_ts", |b| {
        b.iter(|| black_box(gdt_ts(&model, &reference)));
    });
    group.bench_function("gdt_ha", |b| {
        b.iter(|| black_box(gdt_ha(&model, &reference)));
    });
    group.bench_function("weighted_rmsd", |b| {
        let weights = vec![1.0; model.len()];
        b.iter(|| black_box(weighted_rmsd(&model, &reference, &weights)));
    });
    group.finish();

    let guide_count = reference.len().min(128);
    let reference_guides = &reference[..guide_count];
    let model_guides = &model[..guide_count];
    c.bench_function("compare_alignment/ce_128", |b| {
        b.iter(|| {
            black_box(ce_align(
                reference_guides,
                model_guides,
                CeOptions::original(),
            ))
        });
    });

    let first_map: Vec<(u32, u32)> = (0_u32..8_192).map(|index| (index, index + 7)).collect();
    let second_map: Vec<(u32, u32)> = (0_u32..8_192)
        .filter(|index| index % 5 != 0)
        .map(|index| (index, index + 7))
        .collect();
    c.bench_function("compare_contacts/similarity_8192", |b| {
        b.iter(|| black_box(contact_map_similarity(&first_map, &second_map)));
    });
}

fn bench_extended_scores(c: &mut Criterion) {
    let structure = structure(Sample::Medium);
    let reference = coordinates(&structure);
    let model = perturbed(&reference, 0.02);
    let options = DockQOptions {
        contact_distance: 5.0,
        ligand_scale: 8.5,
        interface_scale: 1.5,
    };
    let mut group = c.benchmark_group("compare_extended");
    group.bench_function("dockq/A-B", |b| {
        b.iter(|| black_box(dockq(&structure, &structure, "A", "B", options)));
    });
    group.bench_function("qs_score/A-B", |b| {
        b.iter(|| {
            black_box(qs_score(
                &structure,
                &structure,
                "A",
                "B",
                QsOptions::standard(5.0),
            ))
        });
    });

    let mapping = PointMapping::new(
        (0..reference.len().min(128)).map(|index| PointMatch {
            reference: index,
            model: index,
        }),
        reference.len(),
        model.len(),
    )
    .required("comparison mapping fixture failed");
    let alignment = ComparisonAlignment::NotRequired;
    group.bench_function("interface_rmsd/128", |b| {
        b.iter(|| black_box(interface_rmsd(&reference, &model, &mapping, alignment)));
    });
    group.bench_function("pocket_rmsd/128", |b| {
        b.iter(|| black_box(pocket_rmsd(&reference, &model, &mapping, alignment)));
    });
    group.bench_function("measure_mapping/128", |b| {
        b.iter(|| black_box(measure_mapping(&reference, &model, &mapping, alignment)));
    });

    let areas: Vec<_> = (0_u32..4_096)
        .map(|index| ContactArea {
            first: index,
            second: index + 1,
            area: f64::from(index % 17) + 0.5,
        })
        .collect();
    group.bench_function("cad_score/4096", |b| {
        b.iter(|| black_box(cad_score(&areas, &areas)));
    });

    let component = symmetric_component();
    let ligand_reference = [[0.0, 0.0, 0.0], [-1.2, 0.0, 0.0], [1.2, 0.0, 0.0]];
    let ligand_model = [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [-1.2, 0.0, 0.0]];
    group.bench_function("equivalent_atom_mappings/CO2", |b| {
        b.iter(|| black_box(equivalent_atom_mappings(&component, 4)));
    });
    group.bench_function("ligand_symmetry_rmsd/CO2", |b| {
        b.iter(|| {
            black_box(ligand_symmetry_rmsd(
                &ligand_reference,
                &ligand_model,
                &component,
                4,
            ))
        });
    });
    group.finish();

    let provider = mapping_provider();
    let (reference_structure, target_structure) = mapping_structures();
    c.bench_function("compare_mapping/assign_chains", |b| {
        b.iter(|| {
            black_box(pdbiox_compare::assign_chains(
                &reference_structure,
                &target_structure,
                &provider,
                Namespace::Label,
                pdbiox_seq::Scoring::simple(),
                1.0,
            ))
        });
    });
    c.bench_function("compare_mapping/map_sequence", |b| {
        b.iter(|| {
            black_box(map_sequence_to_structure(
                b"AGS",
                &reference_structure,
                &provider,
                Namespace::Label,
                pdbiox_seq::Scoring::simple(),
            ))
        });
    });
}

fn symmetric_component() -> Component {
    Component {
        id: "CO2".into(),
        name: "carbon dioxide".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: Some("CO2".into()),
        atoms: Arc::from([
            component_atom("C", Element::CARBON),
            component_atom("O1", Element::OXYGEN),
            component_atom("O2", Element::OXYGEN),
        ]),
        bonds: Arc::from([
            ComponentBond {
                atom_a: "C".into(),
                atom_b: "O1".into(),
                order: BondOrder::Double,
                aromatic: false,
                stereo: None,
            },
            ComponentBond {
                atom_a: "C".into(),
                atom_b: "O2".into(),
                order: BondOrder::Double,
                aromatic: false,
                stereo: None,
            },
        ]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn component_atom(name: &str, element: Element) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    }
}

fn mapping_provider() -> MemoryProvider {
    MemoryProvider::new(
        DictionaryVersion::new("compare-bench"),
        [
            component_code("ALA", b'A'),
            component_code("GLY", b'G'),
            component_code("SER", b'S'),
        ],
    )
    .required("comparison provider fixture failed")
}

fn component_code(id: &str, code: u8) -> Component {
    Component {
        id: id.into(),
        name: id.into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(code),
        formula: None,
        atoms: Arc::from([]),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn mapping_structures() -> (pdbiox_core::Structure, pdbiox_core::Structure) {
    const HEADER: &str = "data_bench\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
    let reference = format!(
        "{HEADER}ATOM 1 C CA ALA A 1 0 0 0\nATOM 2 C CA GLY A 2 1 0 0\nATOM 3 C CA SER A 3 2 0 0\n"
    );
    let target = reference.replace(" A ", " B ");
    let read = |text: String| {
        pdbiox_cif::read(
            &pdbiox_bench::input(text.as_bytes()),
            &pdbiox_core::io::ReadOptions::new(),
        )
        .map(|(structure, _)| structure)
        .required("comparison mapping structure failed")
    };
    (read(reference), read(target))
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_scores(&mut criterion);
    bench_extended_scores(&mut criterion);
    criterion.final_summary();
}
