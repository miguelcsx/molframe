use super::{gdt_ha, gdt_ts, gdt_with_cutoffs};

/// A deterministic self-avoiding-ish random walk with 3.8 Å steps, so the
/// fixture looks like a C-alpha trace without depending on a random crate.
fn trace(length: usize) -> Vec<[f32; 3]> {
    let mut state = 0x2545_f491_u32;
    let mut next = move || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        f32::from(u16::try_from(state >> 16).expect("sixteen bits fit")) / 65_535.0 - 0.5
    };
    let mut points = Vec::with_capacity(length);
    let mut at = [0.0f32; 3];
    for _ in 0..length {
        points.push(at);
        let mut step = [next(), next(), next()];
        let norm = (step[0] * step[0] + step[1] * step[1] + step[2] * step[2]).sqrt();
        for axis in &mut step {
            *axis *= 3.8 / norm;
        }
        // Bias the walk forward so it stays extended and the domains separate.
        at = [at[0] + step[0] + 1.5, at[1] + step[1], at[2] + step[2]];
    }
    points
}

/// Rotates a range of points about the z axis through the origin, then shifts
/// them: one rigid motion, applied to one domain.
fn move_domain(points: &mut [[f32; 3]], degrees: f32, shift: [f32; 3]) {
    let (sin, cos) = degrees.to_radians().sin_cos();
    for point in points {
        let [x, y, z] = *point;
        *point = [
            cos * x - sin * y + shift[0],
            sin * x + cos * y + shift[1],
            z + shift[2],
        ];
    }
}

fn plain_fit_gdt(model: &[[f32; 3]], reference: &[[f32; 3]], cutoffs: &[f64]) -> f64 {
    let fit = molframe_geom::superpose(model, reference).expect("valid fit");
    let mut total = 0.0;
    for cutoff in cutoffs {
        let inside = model
            .iter()
            .zip(reference)
            .filter(|&(&point, &target)| {
                molframe_geom::distance(fit.transform.apply(point), target) <= *cutoff
            })
            .count();
        total += f64::from(u32::try_from(inside).expect("small")) / reference.len() as f64;
    }
    total / cutoffs.len() as f64
}

#[test]
fn a_rigidly_displaced_domain_leaves_the_rest_fully_counted() {
    // Sixty untouched atoms plus forty moved by one rigid motion far larger than
    // every cutoff. A superposition that holds the sixty exactly puts them inside
    // every cutoff, and none can hold both groups, so the maximum is exactly 0.6
    // at all four cutoffs. That value follows from the construction alone.
    let reference = trace(100);
    let mut model = reference.clone();
    move_domain(&mut model[60..], 35.0, [30.0, -20.0, 25.0]);

    let ts = gdt_ts(&model, &reference).expect("valid");
    let ha = gdt_ha(&model, &reference).expect("valid");
    assert!((ts - 0.6).abs() < 1e-9, "GDT-TS {ts}");
    assert!((ha - 0.6).abs() < 1e-9, "GDT-HA {ha}");

    let plain = plain_fit_gdt(&model, &reference, &[1.0, 2.0, 4.0, 8.0]);
    assert!(
        plain < 0.5,
        "the single least-squares fit must fall short of the maximum: {plain}"
    );
}

#[test]
fn the_largest_rigid_domain_sets_the_score() {
    // Three domains of 50, 30 and 20 atoms, the last two each moved differently.
    let reference = trace(100);
    let mut model = reference.clone();
    move_domain(&mut model[50..80], -50.0, [25.0, 30.0, -15.0]);
    move_domain(&mut model[80..], 70.0, [-40.0, 10.0, 35.0]);
    let ts = gdt_ts(&model, &reference).expect("valid");
    assert!((ts - 0.5).abs() < 1e-9, "GDT-TS {ts}");
}

#[test]
fn the_search_never_scores_below_the_plain_least_squares_fit() {
    let reference = trace(80);
    let mut model = reference.clone();
    // Spread a deterministic wobble of up to about 1.5 Å over every atom.
    for (index, point) in model.iter_mut().enumerate() {
        let phase = f32::from(u16::try_from(index).expect("small")) * 0.7;
        point[0] += 1.2 * phase.sin();
        point[1] += 0.9 * (phase * 1.3).cos();
        point[2] += 0.6 * (phase * 0.4).sin();
    }
    for cutoffs in [[1.0, 2.0, 4.0, 8.0], [0.5, 1.0, 2.0, 4.0]] {
        let searched = gdt_with_cutoffs(&model, &reference, &cutoffs).expect("valid");
        let plain = plain_fit_gdt(&model, &reference, &cutoffs);
        assert!(searched + 1e-12 >= plain, "{searched} < {plain}");
    }
}

#[test]
fn a_rigid_motion_of_the_whole_model_scores_one() {
    let reference = trace(60);
    let mut model = reference.clone();
    move_domain(&mut model, 120.0, [5.0, -9.0, 14.0]);
    let ts = gdt_ts(&model, &reference).expect("valid");
    assert!((ts - 1.0).abs() < 1e-9, "{ts}");
}
