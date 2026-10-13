use super::{DEFAULT_POINTS, SasaMethod, SasaOptions, SasaWordError, measure_atoms};
use molframe_chem::{HeteroPolicy, MetalFallback, RadiusSet};
use molframe_core::ExecutionContext;

fn areas(radii: &[f32], options: &SasaOptions) -> Vec<f64> {
    let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let Ok(areas) = measure_atoms(&positions, radii, options, &ExecutionContext::default()) else {
        panic!("two atoms are a valid input");
    };
    areas
}

#[test]
fn the_default_is_polymer_only_shrake_rupley_with_bondi_radii() {
    let options = SasaOptions::default();
    assert_eq!(
        options.method,
        SasaMethod::ShrakeRupley {
            points: DEFAULT_POINTS
        }
    );
    assert_eq!(options.hetero, HeteroPolicy::POLYMER);
    assert_eq!(options.metals, MetalFallback::Vdw);
    assert_eq!(options.radii, RadiusSet::Bondi);
    assert!((options.probe - 1.4).abs() < f32::EPSILON);
}

#[test]
fn a_zero_radius_atom_neither_has_an_area_nor_occludes() {
    let options = SasaOptions::default();
    let alone = areas(&[1.7, f32::NAN], &options);
    let beside_zero = areas(&[1.7, 0.0], &options);
    assert!(beside_zero[1].is_nan());
    assert_eq!(alone[0].to_bits(), beside_zero[0].to_bits());
}

#[test]
fn an_infinite_radius_is_left_out_the_same_way() {
    let options = SasaOptions::default();
    let alone = areas(&[1.7, f32::NAN], &options);
    let beside_infinite = areas(&[1.7, f32::INFINITY], &options);
    assert!(beside_infinite[1].is_nan());
    assert_eq!(alone[0].to_bits(), beside_infinite[0].to_bits());
}

#[test]
fn both_methods_measure_the_same_atoms() {
    let lee = SasaOptions {
        method: SasaMethod::LeeRichards { slices: 20 },
        ..SasaOptions::default()
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
    ) else {
        panic!("valid words");
    };
    assert_eq!(options.method, SasaMethod::LeeRichards { slices: 30 });
    assert!(options.hetero.waters && !options.hetero.ligands);
    assert_eq!(options.metals, MetalFallback::Ionic);
    let Ok(defaults) = SasaOptions::from_words(None, None, None, None, None) else {
        panic!("valid words");
    };
    assert_eq!(defaults, SasaOptions::default());
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
    assert!(SasaOptions::from_words(None, None, None, Some("lipids"), None).is_err());
}
