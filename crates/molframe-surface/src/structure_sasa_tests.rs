use super::{
    DEFAULT_POINTS, SasaMethod, SasaOptions, SasaWordError, measure_atoms, structure_sasa_report,
};
use molframe_bench::table_store;
use molframe_chem::{HeteroPolicy, MetalFallback, TableRequest};
use molframe_core::ExecutionContext;

fn bondi_options() -> SasaOptions {
    match SasaOptions::resolve(&TableRequest::named("bondi"), table_store()) {
        Ok(options) => options,
        Err(error) => panic!("bondi radii: {error}"),
    }
}

fn areas(radii: &[f32], options: &SasaOptions) -> Vec<f64> {
    let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let Ok(areas) = measure_atoms(&positions, radii, options, &ExecutionContext::default()) else {
        panic!("two atoms are a valid input");
    };
    areas
}

#[test]
fn the_default_is_polymer_only_shrake_rupley_with_bondi_radii() {
    let options = bondi_options();
    assert_eq!(
        options.method,
        SasaMethod::ShrakeRupley {
            points: DEFAULT_POINTS
        }
    );
    assert_eq!(options.hetero, HeteroPolicy::POLYMER);
    assert_eq!(options.metals, MetalFallback::Vdw);
    assert_eq!(options.radii.identity().name, "bondi");
    assert!((options.probe - 1.4).abs() < f32::EPSILON);
}

#[test]
fn a_zero_radius_atom_neither_has_an_area_nor_occludes() {
    let options = bondi_options();
    let alone = areas(&[1.7, f32::NAN], &options);
    let beside_zero = areas(&[1.7, 0.0], &options);
    assert!(beside_zero[1].is_nan());
    assert_eq!(alone[0].to_bits(), beside_zero[0].to_bits());
}

#[test]
fn an_infinite_radius_is_left_out_the_same_way() {
    let options = bondi_options();
    let alone = areas(&[1.7, f32::NAN], &options);
    let beside_infinite = areas(&[1.7, f32::INFINITY], &options);
    assert!(beside_infinite[1].is_nan());
    assert_eq!(alone[0].to_bits(), beside_infinite[0].to_bits());
}

#[test]
fn both_methods_measure_the_same_atoms() {
    let lee = SasaOptions {
        method: SasaMethod::LeeRichards { slices: 20 },
        ..bondi_options()
    };
    let measured = areas(&[1.7, f32::NAN], &lee);
    assert!(measured[0].is_finite() && measured[1].is_nan());
}

#[test]
fn words_build_the_same_options_as_the_fields() {
    let Ok(options) = SasaOptions::from_words(
        Some("lee-richards"),
        Some(30),
        Some(1.2),
        Some("waters"),
        Some("ionic"),
        &TableRequest::named("bondi"),
        table_store(),
    ) else {
        panic!("valid words");
    };
    assert_eq!(options.method, SasaMethod::LeeRichards { slices: 30 });
    assert!(options.hetero.waters && !options.hetero.ligands);
    assert_eq!(options.metals, MetalFallback::Ionic);
    let Ok(defaults) = SasaOptions::from_words(
        None,
        None,
        None,
        None,
        None,
        &TableRequest::named("bondi"),
        table_store(),
    ) else {
        panic!("valid words");
    };
    assert_eq!(defaults, bondi_options());
}

#[test]
fn bad_words_and_a_zero_count_are_refused() {
    assert!(matches!(
        SasaMethod::named("fast", None),
        Err(SasaWordError::Method(_))
    ));
    assert!(matches!(
        SasaMethod::named("sr", Some(0)),
        Err(SasaWordError::Samples)
    ));
    assert!(
        SasaOptions::from_words(
            None,
            None,
            None,
            Some("lipids"),
            None,
            &TableRequest::named("bondi"),
            table_store()
        )
        .is_err()
    );
}

#[test]
fn a_supplied_table_drives_the_areas_and_is_named_in_the_report() {
    let crambin = molframe_bench::structure(molframe_bench::Sample::Small);
    let positions = molframe_bench::coordinates(&crambin);
    let context = ExecutionContext::default();
    let Ok(default) = structure_sasa_report(&crambin, &positions, &bondi_options(), &context)
    else {
        panic!("default areas");
    };
    assert_eq!(default.policy.radii.name, "bondi");
    let Ok(other) = SasaOptions::resolve(&TableRequest::named("alvarez"), table_store()) else {
        panic!("alvarez radii");
    };
    let Ok(changed) = structure_sasa_report(&crambin, &positions, &other, &context) else {
        panic!("alvarez areas");
    };
    assert_eq!(changed.policy.radii.name, "alvarez");
    assert_ne!(changed.areas, default.areas);
}

#[test]
fn an_unknown_table_names_how_to_supply_it() {
    let refused = SasaOptions::resolve(&TableRequest::named("nowhere"), table_store());
    match refused {
        Err(error) => assert!(error.to_string().contains("nowhere"), "{error}"),
        Ok(_) => panic!("a table nobody supplied was found"),
    }
}
