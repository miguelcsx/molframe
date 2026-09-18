use super::*;

#[test]
fn the_fixed_point_path_matches_the_general_float_parser() {
    let values = [
        "0.0",
        "-0.0",
        "1.0",
        "-1.0",
        "12.345",
        "-12.345",
        "0.001",
        "-0.001",
        "999.999",
        "-999.999",
        "1234.5678",
        "-1234.5678",
        "100.00",
        "0.5",
        "37.50",
        "20.00",
        "1.00",
        "0.000",
        "-0.000",
        "+3.25",
    ];
    for text in values {
        let Some(fast) = fixed_point(text) else {
            panic!("the fixed-point path refused {text}")
        };
        let Ok(general) = text.parse::<f64>() else {
            panic!("the general parser refused {text}")
        };
        assert_eq!(fast.to_bits(), general.to_bits(), "differs for {text}");
    }
}

#[test]
fn a_swept_coordinate_range_matches_the_general_float_parser() {
    let mut value = -500.0_f64;
    while value < 500.0 {
        let text = format!("{value:.3}");
        let Some(fast) = fixed_point(&text) else {
            panic!("the fixed-point path refused {text}")
        };
        let Ok(general) = text.parse::<f64>() else {
            panic!("the general parser refused {text}")
        };
        assert_eq!(fast.to_bits(), general.to_bits(), "differs for {text}");
        value += 0.037;
    }
}

#[test]
fn shapes_outside_the_fast_path_are_declined_rather_than_guessed() {
    for text in [
        "1e5",
        "1E5",
        "1.5e-3",
        "inf",
        "NaN",
        "",
        "-",
        ".",
        "1.2.3",
        "12a.5",
        "0x10",
        "1234567890123456.0",
    ] {
        assert!(
            fixed_point(text).is_none(),
            "the fixed-point path should decline {text}"
        );
    }
}
