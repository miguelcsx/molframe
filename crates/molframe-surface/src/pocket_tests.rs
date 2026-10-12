use super::{PocketScanError, PocketScanOptions, PointProximity, scan_pockets};
use crate::numeric::f64_to_f32;
use crate::sampling::fibonacci_sphere;
use molframe_core::CancellationToken;

/// Atoms on a sphere shell of `radius`, with a round opening towards +z when `opening`.
fn shell(radius: f64, opening: bool) -> (Vec<[f32; 3]>, Vec<f32>) {
    let positions: Vec<[f32; 3]> = fibonacci_sphere(400)
        .into_iter()
        .filter(|d| !opening || d[2] < 0.8)
        .map(|d| [d[0], d[1], d[2]].map(|c| f64_to_f32(c * radius)))
        .collect();
    let radii = vec![1.7; positions.len()];
    (positions, radii)
}

fn scan(
    opening: bool,
    cancel: &CancellationToken,
) -> Result<Vec<super::GridPocket>, PocketScanError> {
    let (positions, radii) = shell(9.0, opening);
    scan_pockets(&positions, &radii, 5, &PocketScanOptions::default(), cancel)
}

#[test]
fn a_shell_with_an_opening_holds_one_pocket_near_its_centre() {
    let Ok(found) = scan(true, &CancellationToken::new()) else {
        panic!("scanned");
    };
    assert!(!found.is_empty());
    let pocket = &found[0];
    assert!(pocket.points.len() >= 24);
    assert!(pocket.centre.iter().all(|c| c.abs() < 4.0), "{:?}", pocket.centre);
    assert!(pocket.volume > 0.0);
}

#[test]
fn a_sealed_shell_holds_no_pocket_because_a_probe_cannot_reach_inside() {
    let Ok(found) = scan(false, &CancellationToken::new()) else {
        panic!("scanned");
    };
    assert!(found.is_empty());
}

#[test]
fn a_cancelled_scan_stops_and_says_so() {
    let token = CancellationToken::new();
    token.cancel();
    assert_eq!(scan(true, &token), Err(PocketScanError::Cancelled));
}

#[test]
fn too_few_atoms_are_refused_by_name() {
    let result = scan_pockets(
        &[[0.0; 3]; 3],
        &[1.7; 3],
        5,
        &PocketScanOptions::default(),
        &CancellationToken::new(),
    );
    assert_eq!(result, Err(PocketScanError::TooFewAtoms { atoms: 3 }));
}

#[test]
fn a_grid_over_the_cell_limit_is_refused_with_its_size() {
    let (positions, radii) = shell(9.0, true);
    let options = PocketScanOptions {
        max_cells: 100,
        ..PocketScanOptions::default()
    };
    let result = scan_pockets(&positions, &radii, 5, &options, &CancellationToken::new());
    assert!(matches!(
        result,
        Err(PocketScanError::TooLarge { limit: 100, .. })
    ));
}

#[test]
fn mismatched_positions_and_radii_are_refused() {
    let result = scan_pockets(
        &[[0.0; 3]; 20],
        &[1.7; 19],
        5,
        &PocketScanOptions::default(),
        &CancellationToken::new(),
    );
    assert_eq!(
        result,
        Err(PocketScanError::LengthMismatch {
            positions: 20,
            radii: 19
        })
    );
}

#[test]
fn a_non_positive_spacing_is_refused() {
    let (positions, radii) = shell(9.0, true);
    let options = PocketScanOptions {
        spacing: 0.0,
        ..PocketScanOptions::default()
    };
    let result = scan_pockets(&positions, &radii, 5, &options, &CancellationToken::new());
    assert_eq!(result, Err(PocketScanError::InvalidOptions));
}

#[test]
fn the_nearest_point_agrees_with_a_brute_force_search_within_the_radius() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        f64_to_f32(f64::from(u32::try_from(state >> 40).unwrap_or(0)) / 16_777_216.0)
    };
    let points: Vec<[f32; 3]> = (0..60)
        .map(|_| [next() * 24.0 - 12.0, next() * 24.0 - 12.0, next() * 24.0 - 12.0])
        .collect();
    let radius = 3.6;
    let index = PointProximity::new(&points, radius);
    for _ in 0..200 {
        let q = [next() * 24.0 - 12.0, next() * 24.0 - 12.0, next() * 24.0 - 12.0];
        let brute = points
            .iter()
            .map(|p| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt())
            .fold(f32::MAX, f32::min);
        let got = index.nearest(q);
        if brute <= radius {
            assert!((got - brute).abs() < f32::EPSILON, "{got} vs {brute}");
        } else {
            assert!(got > radius);
        }
    }
    assert!(!index.may_reach([100.0, 0.0, 0.0], radius));
}
