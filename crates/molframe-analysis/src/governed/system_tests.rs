//! The analysed system: assemblies and crystal contacts reach the kernel.

use crate::GovernedAnalysisError;
use crate::{
    analyse_structure, analyse_trajectory, contacts_kernel, native_contact_fraction_kernel,
};
use molframe_bench::{Sample, structure_from_cif};
use molframe_cif::{parse, read};
use molframe_core::contract::{
    AnalysisPolicy, AssemblyChoice, HydrogenPolicy, PolicyField, SymmetryPolicy,
};
use molframe_core::{Code, ExecutionContext, InputBuffer, ReadOptions, Structure};
use molframe_spatial::SpatialBackend;
use molframe_traj::{Frame, Trajectory};
use molframe_xtal::{ASSEMBLIES_EXTENSION, SYMMETRY_EXTENSION, lower_assemblies, lower_symmetry};

/// One carbon in chain A at the origin; the assembly is that chain and a copy of it moved
/// 3 A along x, so the unit has one atom and the assembly has two.
const DIMER: &str = r"data_dimer
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 1 A 0 0 0
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 0 1 0 0 0 0 1 0
2 1 0 0 3 0 1 0 0 0 0 1 0
_pdbx_struct_assembly.id 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '1,2' 'A'
";

fn dimer() -> Structure {
    let input = InputBuffer::from_bytes(DIMER.as_bytes().to_vec());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    let assemblies = match lower_assemblies(&document) {
        Ok(assemblies) => assemblies,
        Err(findings) => panic!("assembly lowering failed: {findings:?}"),
    };
    let structure = match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("structure lowering failed: {findings:?}"),
    };
    structure.with_extension(ASSEMBLIES_EXTENSION, assemblies)
}

fn crambin() -> Structure {
    let Some(bytes) = Sample::Tiny.cif() else {
        panic!("1CRN ships an mmCIF");
    };
    let input = InputBuffer::from_bytes(bytes.to_vec());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    let symmetry = match lower_symmetry(&document) {
        Ok(symmetry) => symmetry,
        Err(findings) => panic!("symmetry failed: {findings:?}"),
    };
    let Some(structure) = structure_from_cif(Sample::Tiny) else {
        panic!("1CRN reads");
    };
    structure.with_extension(SYMMETRY_EXTENSION, symmetry)
}

fn policy(assembly: AssemblyChoice, symmetry: SymmetryPolicy) -> AnalysisPolicy {
    AnalysisPolicy {
        assembly,
        symmetry,
        ..AnalysisPolicy::default()
    }
}

fn contact_count(structure: &Structure, policy: &AnalysisPolicy, cutoff: f32) -> usize {
    let kernel = contacts_kernel(cutoff, SpatialBackend::BruteForce);
    match analyse_structure(structure, policy, &kernel, &ExecutionContext::default()) {
        Ok(result) => result.value().map_or(0, crate::ContactTable::len),
        Err(error) => panic!("analysis failed: {error}"),
    }
}

#[test]
fn the_unit_has_no_contact_with_itself_and_its_assembly_does() {
    let structure = dimer();
    let unit = AnalysisPolicy::default();
    let assembly = policy(AssemblyChoice::Biological("1".into()), SymmetryPolicy::None);
    assert_eq!(contact_count(&structure, &unit, 4.0), 0);
    assert_eq!(contact_count(&structure, &assembly, 4.0), 1);
    assert_eq!(contact_count(&structure, &assembly, 2.0), 0);
}

#[test]
fn the_assembly_records_what_it_is_and_what_it_read() {
    let structure = dimer();
    let assembly = policy(
        AssemblyChoice::Biological("1".into()),
        SymmetryPolicy::BiologicalAssembly,
    );
    let kernel = contacts_kernel(4.0, SpatialBackend::BruteForce);
    let Ok(result) =
        analyse_structure(&structure, &assembly, &kernel, &ExecutionContext::default())
    else {
        panic!("the assembly analysis should run");
    };
    assert_eq!(result.coverage.intended, 2);
    assert!(
        result
            .assumptions
            .iter()
            .any(|assumption| assumption.field == PolicyField::Assembly
                && assumption.value.contains("biological assembly 1"))
    );
    let Some(reads) = result.provenance.policy_reads() else {
        panic!("the executor records what it applied");
    };
    for field in [
        PolicyField::Assembly,
        PolicyField::Altloc,
        PolicyField::Hydrogens,
    ] {
        assert!(reads.contains(&field));
    }
    assert!(
        !reads.contains(&PolicyField::VdwRadii),
        "contacts never read radii"
    );
}

#[test]
fn a_trajectory_of_the_unit_moves_its_copies_with_it() {
    let structure = dimer();
    let assembly = policy(AssemblyChoice::Biological("1".into()), SymmetryPolicy::None);
    let frames = vec![
        Frame {
            positions: vec![[0.0, 0.0, 0.0]],
        },
        Frame {
            positions: vec![[0.0, 5.0, 0.0]],
        },
    ];
    let Ok(trajectory) = Trajectory::from_frames(frames) else {
        panic!("two one-atom frames");
    };
    let kernel = contacts_kernel(4.0, SpatialBackend::BruteForce);
    let Ok(result) = analyse_trajectory(
        &structure,
        &trajectory,
        &assembly,
        &kernel,
        &ExecutionContext::default(),
    ) else {
        panic!("trajectory should run");
    };
    // The copy sits 3 A from the unit in every frame, because it is the unit moved by the
    // same rigid motion; it is not a second set of coordinates.
    let distances: Vec<Vec<u32>> = result
        .value()
        .map(|tables| {
            tables
                .iter()
                .map(|table| table.distances().iter().map(|d| d.to_bits()).collect())
                .collect()
        })
        .unwrap_or_default();
    assert_eq!(
        distances,
        vec![vec![3.0_f32.to_bits()], vec![3.0_f32.to_bits()]]
    );
}

#[test]
fn a_system_that_cannot_be_built_says_why_before_anything_runs() {
    let kernel = contacts_kernel(4.0, SpatialBackend::BruteForce);
    let context = ExecutionContext::default();
    let dimer = dimer();
    let Err(missing) = analyse_structure(
        &dimer,
        &policy(AssemblyChoice::Biological("9".into()), SymmetryPolicy::None),
        &kernel,
        &context,
    ) else {
        panic!("an assembly the file lacks must be refused");
    };
    assert!(matches!(&missing, GovernedAnalysisError::System(f) if f[0].code() == Code::E6002));

    let Err(contradiction) = analyse_structure(
        &dimer,
        &policy(
            AssemblyChoice::AsymmetricUnit,
            SymmetryPolicy::Crystallographic,
        ),
        &kernel,
        &context,
    ) else {
        panic!("a contradictory pair must be refused");
    };
    assert!(
        matches!(&contradiction, GovernedAnalysisError::System(f) if f[0].code() == Code::E6004)
    );

    let Err(no_cell) = analyse_structure(
        &dimer,
        &policy(
            AssemblyChoice::Crystal { radius: 4.0 },
            SymmetryPolicy::None,
        ),
        &kernel,
        &context,
    ) else {
        panic!("crystal contacts need symmetry operators");
    };
    assert!(matches!(&no_cell, GovernedAnalysisError::System(f) if f[0].code() == Code::E6016));
}

#[test]
fn an_analysis_that_pairs_atoms_across_structures_refuses_copies() {
    let structure = dimer();
    let target = dimer();
    let kernel = native_contact_fraction_kernel(&target, 4.0, 1.2, SpatialBackend::BruteForce);
    let assembly = policy(AssemblyChoice::Biological("1".into()), SymmetryPolicy::None);
    let refused = analyse_structure(&structure, &assembly, &kernel, &ExecutionContext::default());
    assert!(matches!(
        refused,
        Err(GovernedAnalysisError::ReplicatedSystemUnsupported(_))
    ));
}

#[test]
fn hydrogens_the_policy_excludes_are_not_analysed_and_inference_is_refused() {
    let text = DIMER.replace(
        "1 C CA GLY A 1 1 A 0 0 0\n",
        "1 C CA GLY A 1 1 A 0 0 0\n2 H H1 GLY A 1 1 A 1 0 0\n",
    );
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let structure = match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("structure lowering failed: {findings:?}"),
    };
    let explicit = AnalysisPolicy::default();
    let excluded = AnalysisPolicy {
        hydrogens: HydrogenPolicy::Exclude,
        ..AnalysisPolicy::default()
    };
    assert_eq!(contact_count(&structure, &explicit, 2.0), 1);
    assert_eq!(contact_count(&structure, &excluded, 2.0), 0);

    let inferred = AnalysisPolicy {
        hydrogens: HydrogenPolicy::IncludeInferred,
        ..AnalysisPolicy::default()
    };
    let kernel = contacts_kernel(2.0, SpatialBackend::BruteForce);
    assert!(matches!(
        analyse_structure(&structure, &inferred, &kernel, &ExecutionContext::default()),
        Err(GovernedAnalysisError::UnsupportedPolicyValue("hydrogens"))
    ));
}

#[test]
fn crystal_contacts_of_crambin_equal_the_pairs_a_brute_force_finds_between_unit_and_mates() {
    let structure = crambin();
    let (radius, cutoff) = (4.0_f32, 3.5_f32);
    let unit = structure.atom_count() as usize;
    let crystal = policy(
        AssemblyChoice::Crystal { radius },
        SymmetryPolicy::Crystallographic,
    );
    let kernel = contacts_kernel(cutoff, SpatialBackend::BruteForce);
    let Ok(result) = analyse_structure(&structure, &crystal, &kernel, &ExecutionContext::default())
    else {
        panic!("crambin has crystal contacts");
    };
    let Some(table) = result.value() else {
        panic!("a determinate result");
    };
    let cross = table
        .first()
        .iter()
        .zip(table.second())
        .filter(|(a, b)| (a.as_usize() < unit) != (b.as_usize() < unit))
        .count();
    assert!(
        cross > 0,
        "crambin touches its symmetry mates within {cutoff} A"
    );

    // Independent count: every symmetry operation and lattice translation, whether or not
    // the system kept it, and every pair of a deposited atom with an atom of that image.
    let Some(cell) = structure.data().cell else {
        panic!("1CRN has a cell");
    };
    let Ok(transform) = molframe_xtal::CellTransform::new(&cell) else {
        panic!("a valid cell");
    };
    let Some(symmetry) = molframe_xtal::SymmetryExt::symmetry_set(&structure) else {
        panic!("1CRN carries operators");
    };
    let positions = structure.positions();
    let mut expected = 0;
    for operation in symmetry.operations() {
        for a in -2..=2 {
            for b in -2..=2 {
                for c in -2..=2 {
                    if operation.is_identity() && [a, b, c] == [0, 0, 0] {
                        continue;
                    }
                    let Ok(motion) =
                        molframe_xtal::operation_motion(&transform, operation, [a, b, c])
                    else {
                        panic!("proper operation");
                    };
                    for &image in positions {
                        let moved = motion.apply(image);
                        expected += positions
                            .iter()
                            .filter(|atom| {
                                let d: f32 =
                                    (0..3).map(|axis| (moved[axis] - atom[axis]).powi(2)).sum();
                                d.sqrt() <= cutoff
                            })
                            .count();
                    }
                }
            }
        }
    }
    assert_eq!(
        cross, expected,
        "unit-to-mate pairs, by brute force over every image"
    );

    // The unit alone has none of these: the same cutoff over the deposited atoms.
    let alone = contact_count(&structure, &AnalysisPolicy::default(), cutoff);
    let within_unit = table
        .first()
        .iter()
        .zip(table.second())
        .filter(|(a, b)| a.as_usize() < unit && b.as_usize() < unit)
        .count();
    assert_eq!(
        within_unit, alone,
        "the unit's own contacts are unchanged by its neighbours"
    );
}

#[test]
fn the_result_says_which_input_atom_each_analysed_atom_is() {
    let kernel = contacts_kernel(4.0, SpatialBackend::BruteForce);
    let context = ExecutionContext::default();
    let origin = |structure: &Structure, policy: &AnalysisPolicy| -> Option<Vec<u32>> {
        match analyse_structure(structure, policy, &kernel, &context) {
            Ok(result) => result.atom_origin().map(<[u32]>::to_vec),
            Err(error) => panic!("analysis failed: {error}"),
        }
    };
    // An assembly of one carbon and its copy: two analysed atoms, one input atom.
    let assembly = policy(AssemblyChoice::Biological("1".into()), SymmetryPolicy::None);
    assert_eq!(origin(&dimer(), &assembly), Some(vec![0, 0]));

    // Excluding hydrogens drops the atom that was there: the carbon is input atom 0.
    let text = DIMER.replace(
        "1 C CA GLY A 1 1 A 0 0 0\n",
        "1 H H1 GLY A 1 1 A 1 0 0\n2 C CA GLY A 1 1 A 0 0 0\n",
    );
    let input = InputBuffer::from_bytes(text.into_bytes());
    let Ok((hydrogen_first, _)) = read(&input, &ReadOptions::new()) else {
        panic!("a valid fixture");
    };
    let excluded = AnalysisPolicy {
        hydrogens: HydrogenPolicy::Exclude,
        ..AnalysisPolicy::default()
    };
    assert_eq!(origin(&hydrogen_first, &excluded), Some(vec![1]));
    assert_eq!(
        origin(&hydrogen_first, &AnalysisPolicy::default()),
        Some(vec![0, 1])
    );
}

#[test]
fn the_ceiling_on_candidate_images_is_the_contexts_and_exceeding_it_is_an_error() {
    let structure = crambin();
    let crystal = policy(
        AssemblyChoice::Crystal { radius: 4.0 },
        SymmetryPolicy::Crystallographic,
    );
    let kernel = contacts_kernel(3.5, SpatialBackend::BruteForce);
    let Ok(tight) = ExecutionContext::builder().image_search_limit(10).build() else {
        panic!("a valid context");
    };
    let Err(refused) = analyse_structure(&structure, &crystal, &kernel, &tight) else {
        panic!("ten candidate images cannot cover crambin's neighbourhood");
    };
    assert!(matches!(&refused, GovernedAnalysisError::System(f) if f[0].code() == Code::E6017));
    let Ok(generous) = ExecutionContext::builder()
        .image_search_limit(50_000_000)
        .build()
    else {
        panic!("a valid context");
    };
    assert!(analyse_structure(&structure, &crystal, &kernel, &generous).is_ok());
}
