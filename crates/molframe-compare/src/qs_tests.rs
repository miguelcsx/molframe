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

#[test]
fn a_shifted_contact_matches_the_closed_form() {
    // Two contacts; one moves from 5 Å (weight 1) to 9.28 Å (weight e^-2).
    let native = complex(5.0, 5.0);
    let model = complex(5.0 + 4.28, 5.0);
    let w = (-2.0_f64).exp();
    let expected = (w + 1.0) / (0.5 * ((1.0 + 1.0) + (w * w + 1.0)));
    assert!((qs(&model, &native) - expected).abs() < 1e-6);
}

#[test]
fn an_extra_model_contact_lowers_the_score_by_its_weight() {
    // Native has one contact; the model adds a second of full weight.
    let native = complex(4.0, 40.0);
    let model = complex(4.0, 4.0);
    assert!((qs(&model, &native) - 1.0 / 1.5).abs() < 1e-12);
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
