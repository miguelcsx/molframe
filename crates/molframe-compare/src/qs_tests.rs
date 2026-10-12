use super::{EmptyQsPolicy, QsOptions, contact_weight, qs_score};
use crate::CompareError;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(body: &str) -> Structure {
    let input = InputBuffer::from_bytes(format!("{HEADER}{body}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

/// Chain A residues sit at the origin and 50 Å away; chain B residues sit above
/// them at the given distances. Residues are ALA with a CB, so the CB is the
/// representative atom, placed on the z axis above the origin of its residue.
fn complex(first_gap: f64, second_gap: f64) -> Structure {
    structure(&format!(
        "ATOM 1 C CB ALA A 1 0 0 0\n\
ATOM 2 C CB ALA A 2 50 0 0\n\
ATOM 3 C CB ALA B 1 0 0 {first_gap}\n\
ATOM 4 C CB ALA B 2 50 0 {second_gap}\n"
    ))
}

fn qs(model: &Structure, native: &Structure) -> f64 {
    match qs_score(model, native, "A", "B", QsOptions::published()) {
        Ok(value) => value,
        Err(error) => panic!("valid QS: {error}"),
    }
}

#[test]
fn the_weight_is_flat_then_gaussian() {
    assert!((contact_weight(3.0) - 1.0).abs() < 1e-15);
    assert!((contact_weight(5.0) - 1.0).abs() < 1e-15);
    // One width beyond the plateau: exp(-2).
    assert!((contact_weight(5.0 + 4.28) - (-2.0_f64).exp()).abs() < 1e-12);
}

#[test]
fn an_identical_interface_scores_one() {
    let native = complex(4.0, 9.0);
    assert!((qs(&native, &native) - 1.0).abs() < 1e-12);
}

#[test]
fn an_interface_with_no_shared_contact_scores_zero() {
    let native = complex(4.0, 40.0);
    // The native contact is lost; the model forms a different one.
    let model = structure(
        "ATOM 1 C CB ALA A 1 0 0 0\n\
ATOM 2 C CB ALA A 2 50 0 0\n\
ATOM 3 C CB ALA B 1 0 0 40\n\
ATOM 4 C CB ALA B 2 50 0 4\n",
    );
    assert!(qs(&model, &native).abs() < 1e-12);
}

#[test]
fn a_contact_pushed_further_away_lowers_the_score_monotonically() {
    let native = complex(4.0, 40.0);
    let mut previous = f64::INFINITY;
    for gap in [4.0, 6.0, 8.0, 10.0, 11.5] {
        let value = qs(&complex(gap, 40.0), &native);
        assert!(value < previous, "gap {gap}: {value} !< {previous}");
        previous = value;
    }
    assert!(previous > 0.0);
}

/// `w(min(d1, d2)) * (1 - |d1 - d2| / 12)` for one shared contact, which is the
/// published contribution of that pair to the numerator.
fn shared_term(d_reference: f64, d_model: f64) -> f64 {
    let smaller = if d_reference < d_model {
        d_reference
    } else {
        d_model
    };
    contact_weight(smaller) * (1.0 - (d_reference - d_model).abs() / 12.0)
}

#[test]
fn a_shared_contact_moved_from_four_to_ten_angstrom_scores_one_half() {
    // One contact in both, moved 6 Å. Its term is w(4) * (1 - 6/12) = 0.5 and
    // nothing else is in contact, so the score is exactly the penalty.
    let native = complex(4.0, 40.0);
    let model = complex(10.0, 40.0);
    assert!(
        (qs(&model, &native) - 0.5).abs() < 1e-9,
        "{}",
        qs(&model, &native)
    );
}

#[test]
fn a_shifted_contact_matches_the_published_closed_form() {
    // Two shared contacts, one held at 5 Å and one moved to 5 Å + 4.28 Å.
    // w(5) = 1, so the moved pair contributes 1 - 4.28/12 = 0.643333…;
    // the shared weight sum is 2 and there are no non-shared contacts.
    let native = complex(5.0, 5.0);
    let model = complex(5.0 + 4.28, 5.0);
    // One shared contact is held exactly on its reference and the other is moved
    // 4.28 Å, so the numerator is their two per-pair terms and the denominator is
    // the shared weight of both.
    let numerator = shared_term(5.0, 5.0) + shared_term(5.0, 5.0 + 4.28);
    let expected = numerator / 2.0;
    assert!(
        (expected - 0.821_666_666_666_666_7).abs() < 1e-12,
        "closed form drifted: {expected}"
    );
    let value = qs(&model, &native);
    assert!((value - expected).abs() < 1e-6, "{value} != {expected}");
}

#[test]
fn an_extra_model_contact_halves_the_score_and_swapping_sides_keeps_it() {
    // Native has one contact; the model adds a second of full weight. The
    // non-shared contact enters the denominator alone, so the score is
    // 1 / (1 + 1) rather than a ratio of squared-weight masses.
    let native = complex(4.0, 40.0);
    let model = complex(4.0, 4.0);
    assert!((qs(&model, &native) - 0.5).abs() < 1e-12);
    assert!((qs(&native, &model) - 0.5).abs() < 1e-12);
}

#[test]
fn a_contact_exactly_at_the_boundary_is_part_of_the_interface() {
    // The published contact test is `d <= 12`. One pair at 12 Å on the
    // reference side against 11 Å on the model side is in contact on both, and
    // its weight cancels: the score is the distance penalty alone, 11/12.
    let native = complex(12.0, 40.0);
    let model = complex(11.0, 40.0);
    let value = qs(&model, &native);
    assert!(
        (value - 11.0 / 12.0).abs() < 1e-9,
        "boundary contact dropped: {value}"
    );
}

#[test]
fn glycine_uses_ca_and_other_residues_without_cb_take_no_part() {
    let native = structure(
        "ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 C CA GLY B 1 0 0 4\n",
    );
    assert!((qs(&native, &native) - 1.0).abs() < 1e-12);
    let alanine_without_cb = structure(
        "ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA ALA B 1 0 0 4\n",
    );
    assert!(matches!(
        qs_score(
            &alanine_without_cb,
            &alanine_without_cb,
            "A",
            "B",
            QsOptions::published()
        ),
        Err(CompareError::NoComparablePairs)
    ));
}

#[test]
fn empty_interface_and_cutoff_policies_are_explicit() {
    let native = complex(4.0, 40.0);
    // Chains far apart: representatives exist but nothing is in contact.
    let apart = structure(
        "ATOM 1 C CB ALA A 1 0 0 0\n\
ATOM 2 C CB ALA B 1 0 0 90\n",
    );
    assert!(
        (qs_score(&apart, &apart, "A", "B", QsOptions::published()).unwrap_or(f64::NAN) - 1.0)
            .abs()
            < 1e-12
    );
    assert!(matches!(
        qs_score(
            &apart,
            &apart,
            "A",
            "B",
            QsOptions {
                contact_distance: 12.0,
                empty_policy: EmptyQsPolicy::Error,
            }
        ),
        Err(CompareError::NoComparablePairs)
    ));
    assert!(matches!(
        qs_score(&native, &native, "A", "B", QsOptions::standard(f32::NAN)),
        Err(CompareError::InvalidDistanceCutoff)
    ));
}
