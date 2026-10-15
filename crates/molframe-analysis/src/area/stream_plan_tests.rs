use super::{SasaStreamError, SasaStreamPlan};
use crate::cif_fixture::read_cif;
use molframe_chem::{RadiusSelection, TableRequest, TableStore, VdwRadii};
use molframe_core::ExecutionContext;
use molframe_core::contract::{AnalysisPolicy, PolicyField, RadiiSet};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::Structure;
use molframe_surface::{SasaOptions, measure_atoms};
use molframe_traj::{Frame, StreamingReader};

const TEXT: &str = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N ALA A 1 0 0 0\nATOM 2 C CA ALA A 1 1.5 0 0\n\
ATOM 3 H H ALA A 1 0.5 1 0\nATOM 4 C C GLY A 2 3.0 0 0\n\
HETATM 5 O O HOH B 3 2.0 2 0\n";

fn topology() -> Structure {
    read_cif(TEXT)
}

fn store() -> &'static TableStore {
    molframe_bench::table_store()
}

fn options(table: &str, samples: u16) -> SasaOptions {
    match SasaOptions::from_words(
        None,
        Some(samples),
        None,
        None,
        None,
        &TableRequest::<VdwRadii>::named(table),
        store(),
    ) {
        Ok(options) => options,
        Err(error) => panic!("options: {error}"),
    }
}

fn policy(set: RadiiSet) -> AnalysisPolicy {
    AnalysisPolicy {
        vdw_radii: set,
        ..AnalysisPolicy::default()
    }
}

fn plan(table: &str, policy: &AnalysisPolicy) -> Result<SasaStreamPlan, SasaStreamError> {
    let selection = RadiusSelection::choose(Some(table), true);
    SasaStreamPlan::new(&topology(), &options(table, 64), &selection, policy)
}

#[test]
fn the_plan_measures_the_polymer_heavy_atoms_only() {
    let plan = plan("bondi", &policy(RadiiSet::Bondi)).unwrap();
    assert_eq!(plan.atoms(), &[0, 1, 3]);
}

#[test]
fn streamed_totals_equal_the_batch_areas_summed_over_the_same_atoms() {
    let topology = topology();
    let options = options("bondi", 64);
    let policy = policy(RadiiSet::Bondi);
    let plan = plan("bondi", &policy).unwrap();
    let context = ExecutionContext::default();
    let radii = crate::sasa_atoms::measured_radii(&topology, &options).radii;
    let frames: Vec<Vec<[f32; 3]>> = (0..4u8)
        .map(|step| {
            topology
                .positions()
                .iter()
                .map(|p| [p[0] + 0.1 * f32::from(step), p[1], p[2]])
                .collect()
        })
        .collect();
    let mut reader = StreamingReader::new(
        frames.iter().map(|frame| Frame {
            positions: plan.atoms().iter().map(|&atom| frame[atom]).collect(),
        }),
        plan.atoms().len(),
    );
    let result = plan.run(&mut reader, 5, &context, 4096).unwrap();
    let Some(series) = result.value() else {
        panic!("a determinate series");
    };
    assert_eq!(series.frames, [0, 1, 2, 3]);
    for (frame, total) in frames.iter().zip(&series.areas) {
        let atoms = measure_atoms(frame, &radii, &options, &context).unwrap();
        let batch: f64 = atoms.iter().filter(|area| area.is_finite()).sum();
        assert!((batch - total).abs() < 1e-9, "{batch} vs {total}");
    }
    let applied = result
        .assumptions
        .iter()
        .find(|assumption| assumption.field == PolicyField::VdwRadii);
    assert!(applied.is_some_and(|a| a.value.to_string().contains("bondi")));
}

#[test]
fn a_table_that_conflicts_with_the_policy_is_a_typed_error() {
    let refused = plan("bondi", &policy(RadiiSet::Charmm));
    let Err(SasaStreamError::Radii(error)) = refused else {
        panic!("the conflict is refused");
    };
    assert_eq!(Diagnostic::from(error).code(), Code::E7103);
}

#[test]
fn a_trajectory_of_another_size_is_refused() {
    let plan = plan("bondi", &policy(RadiiSet::Bondi)).unwrap();
    let mut reader = StreamingReader::new(std::iter::empty::<Frame>(), 3);
    let refused = plan.run(&mut reader, 7, &ExecutionContext::default(), 4096);
    assert!(matches!(
        refused,
        Err(SasaStreamError::TopologyAtoms {
            topology: 5,
            trajectory: 7
        })
    ));
}
