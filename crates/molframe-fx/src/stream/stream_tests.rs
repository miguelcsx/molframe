use crate::{
    AtomSite, ComponentRole, ComponentSpec, Constraint, MappingError, MappingPlan,
    MeasurementOptions, Motif, NamedConstraint, align_intrinsic, map_motif, measure_constraints,
};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::execution::CancellationToken;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use std::fmt::Write as _;
use std::sync::Arc;

/// A deterministic scatter of SER, TRP and MET residues in a 30 angstrom box.
fn structure(residues: u32) -> Structure {
    seeded(residues, 12_345)
}

/// The same scatter as [`structure`], placed by `seed`.
fn seeded(residues: u32, seed: u64) -> Structure {
    let mut source = String::from(
        "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f64::from(u32::try_from(state >> 40).unwrap_or(0)) / 16_777_216.0 * 30.0
    };
    let mut serial = 0;
    for residue in 1..=residues {
        let (name, atoms): (&str, &[(&str, &str)]) = match residue % 3 {
            0 => ("SER", &[("OG", "O"), ("CB", "C")]),
            1 => ("TRP", &[("NE1", "N"), ("CA", "C")]),
            _ => ("MET", &[("SD", "S"), ("CA", "C")]),
        };
        let (x, y, z) = (next(), next(), next());
        for (offset, (atom, element)) in atoms.iter().enumerate() {
            serial += 1;
            let shift = f64::from(u32::try_from(offset).unwrap_or(0)) * 1.4;
            let _ = writeln!(
                source,
                "ATOM {serial} {element} {atom} {name} A {residue} {:.3} {y:.3} {z:.3}",
                x + shift
            );
        }
    }
    let input = InputBuffer::from_bytes(source.into_bytes());
    molframe_cif::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("fixture failed: {findings:?}"))
        .0
}

fn residue_motif(components: &[(&str, &str, &str)], distances: &[(&str, &str, f64, f64)]) -> Motif {
    let specs = components.iter().map(|(name, id, atom)| {
        (
            Box::<str>::from(*name),
            ComponentSpec::new(ComponentRole::Residue)
                .component(*id)
                .require(*atom),
        )
    });
    let atom_of = |component: &str| {
        components
            .iter()
            .find(|(name, ..)| *name == component)
            .map_or("", |(.., atom)| *atom)
    };
    let constraints = distances
        .iter()
        .enumerate()
        .map(|(index, (a, b, t, tol))| NamedConstraint {
            name: format!("d{index}").into(),
            constraint: Constraint::Distance {
                first: AtomSite::new(*a, atom_of(a)),
                second: AtomSite::new(*b, atom_of(b)),
                target: *t,
                tolerance: *tol,
            },
        });
    Motif::new(specs, constraints).unwrap_or_else(|error| panic!("motif failed: {error}"))
}

fn stream_all(plan: &Arc<MappingPlan>) -> Vec<Vec<u32>> {
    let cancel = CancellationToken::new();
    let mut stream = plan.stream();
    let mut tuple = vec![0; plan.width()];
    let mut all = Vec::new();
    while stream
        .next_into(&mut tuple, &cancel)
        .unwrap_or_else(|error| panic!("stream failed: {error}"))
    {
        all.push(tuple.clone());
    }
    all
}

fn plan_for(structure: &Structure, motif: &Motif, constrained: bool) -> Arc<MappingPlan> {
    Arc::new(
        MappingPlan::build(
            structure,
            motif,
            None,
            &AnalysisPolicy::default(),
            constrained,
        )
        .unwrap_or_else(|error| panic!("plan failed: {error}")),
    )
}

/// Legacy exhaustive mappings, kept when every distance constraint measures satisfied.
fn reference(structure: &Structure, motif: &Motif, filter: bool) -> Vec<Vec<u32>> {
    let set = map_motif(
        structure,
        motif,
        None,
        &AnalysisPolicy::default(),
        usize::MAX,
    )
    .unwrap_or_else(|error| panic!("mapping failed: {error}"));
    let options = MeasurementOptions {
        maximum_alternatives: 256,
        plane_fit: molframe_geom::EigenOptions::standard(),
    };
    set.mappings
        .into_iter()
        .filter(|mapping| {
            if !filter {
                return true;
            }
            let aligned = align_intrinsic(mapping.clone());
            let measured = measure_constraints(structure, motif, &aligned, options)
                .unwrap_or_else(|error| panic!("measure failed: {error}"));
            measured
                .constraints
                .iter()
                .all(|constraint| constraint.satisfied == Some(true))
        })
        .map(|mapping| mapping.components.values().map(|r| r.get()).collect())
        .collect()
}

#[test]
fn unconstrained_stream_matches_the_exhaustive_order_and_count() {
    let structure = structure(60);
    let motif = residue_motif(
        &[("a", "SER", "OG"), ("b", "TRP", "NE1")],
        &[("a", "b", 6.0, 2.0)],
    );
    let plan = plan_for(&structure, &motif, false);
    let streamed = stream_all(&plan);
    assert_eq!(streamed, reference(&structure, &motif, false));
    assert_eq!(plan.count(&CancellationToken::new()).ok(), Some(400));
}

#[test]
fn constrained_pair_equals_brute_force_filter_in_the_same_order() {
    let structure = structure(240);
    let motif = residue_motif(
        &[("a", "SER", "OG"), ("b", "TRP", "NE1")],
        &[("a", "b", 6.0, 2.0)],
    );
    let plan = plan_for(&structure, &motif, true);
    let expected = reference(&structure, &motif, true);
    assert!(!expected.is_empty());
    assert!(expected.len() < 80 * 80);
    assert_eq!(stream_all(&plan), expected);
    assert_eq!(
        plan.count(&CancellationToken::new()).ok(),
        Some(expected.len() as u64)
    );
}

#[test]
fn constrained_ternary_chain_equals_brute_force() {
    let structure = structure(150);
    let motif = residue_motif(
        &[("a", "SER", "OG"), ("b", "TRP", "NE1"), ("c", "MET", "SD")],
        &[
            ("a", "b", 7.0, 3.0),
            ("b", "c", 8.0, 3.0),
            ("c", "a", 9.0, 4.0),
        ],
    );
    let plan = plan_for(&structure, &motif, true);
    let expected = reference(&structure, &motif, true);
    assert!(!expected.is_empty());
    assert_eq!(stream_all(&plan), expected);
    assert_eq!(
        plan.count(&CancellationToken::new()).ok(),
        Some(expected.len() as u64)
    );
}

#[test]
fn a_residue_never_fills_two_components() {
    let structure = structure(30);
    let motif = residue_motif(
        &[("a", "SER", "OG"), ("b", "SER", "OG")],
        &[("a", "b", 10.0, 10.0)],
    );
    let plan = plan_for(&structure, &motif, true);
    let streamed = stream_all(&plan);
    assert!(streamed.iter().all(|tuple| tuple[0] != tuple[1]));
    assert_eq!(streamed, reference(&structure, &motif, true));
    let free = plan_for(&structure, &motif, false);
    assert_eq!(free.count(&CancellationToken::new()).ok(), Some(90));
}

#[test]
fn a_cancelled_token_stops_enumeration() {
    let structure = structure(9000);
    let motif = residue_motif(&[("a", "SER", "OG"), ("b", "TRP", "NE1")], &[]);
    let plan = plan_for(&structure, &motif, false);
    let cancel = CancellationToken::new();
    cancel.cancel();
    assert!(matches!(plan.count(&cancel), Err(MappingError::Cancelled)));
    let mut stream = plan.stream();
    let mut tuple = vec![0; 2];
    let mut outcome = Ok(true);
    while matches!(outcome, Ok(true)) {
        outcome = stream.next_into(&mut tuple, &cancel);
    }
    assert!(matches!(outcome, Err(MappingError::Cancelled)));
}

#[test]
fn the_limit_message_speaks_of_mappings() {
    let structure = structure(60);
    let motif = residue_motif(&[("a", "SER", "OG"), ("b", "TRP", "NE1")], &[]);
    let error = map_motif(&structure, &motif, None, &AnalysisPolicy::default(), 5)
        .err()
        .map(|error| error.to_string())
        .unwrap_or_default();
    assert!(error.contains("limit of 5 mappings"), "{error}");
}

#[path = "property_tests.rs"]
mod property;
