use super::*;

#[test]
fn rigidly_moved_chain_aligns_every_guide_atom() {
    let reference: Vec<[f32; 3]> = (0_i16..24)
        .map(|index| {
            let value = f32::from(index);
            [value, (value * 0.4).sin() * 3.0, (value * 0.7).cos() * 2.0]
        })
        .collect();
    let mobile: Vec<[f32; 3]> = reference
        .iter()
        .map(|point| [-point[1] + 10.0, point[0] - 4.0, point[2] + 2.0])
        .collect();
    let alignment = ce_align(&reference, &mobile, CeOptions::original())
        .unwrap_or_else(|error| panic!("CE failed: {error}"));
    assert_eq!(alignment.reference_indices.len(), 24);
    assert_eq!(alignment.reference_indices, alignment.mobile_indices);
    assert!(alignment.rmsd < 1.0e-5);
}

#[test]
fn insertion_is_crossed_by_a_bounded_gap() {
    let reference: Vec<[f32; 3]> = (0_i16..24)
        .map(|index| {
            [
                f32::from(index),
                (f32::from(index * index) * 0.03).sin(),
                f32::from(index % 5),
            ]
        })
        .collect();
    let mut mobile = reference[..8].to_vec();
    mobile.push([100.0, 100.0, 100.0]);
    mobile.extend_from_slice(&reference[8..]);
    let alignment = ce_align(
        &reference,
        &mobile,
        CeOptions {
            max_gap: 2,
            ..CeOptions::original()
        },
    )
    .unwrap_or_else(|error| panic!("gapped CE failed: {error}"));
    assert!(alignment.reference_indices.len() >= 16);
    assert!(alignment.mobile_indices.iter().any(|index| *index >= 9));
}

#[test]
fn two_fragments_are_required() {
    let points = vec![[0.0; 3]; 15];
    assert!(matches!(
        ce_align(&points, &points, CeOptions::original()),
        Err(CeError::TooFewPoints {
            required: 16,
            actual: 15
        })
    ));
}

#[test]
fn thresholds_and_significance_are_explicit_and_validated() {
    let points: Vec<_> = (0_i16..24)
        .map(|index| [f32::from(index), (f32::from(index) * 0.3).sin(), 0.0])
        .collect();
    let alignment = ce_align(
        &points,
        &points,
        CeOptions {
            significance: None,
            ..CeOptions::original()
        },
    )
    .unwrap_or_else(|error| panic!("explicit CE failed: {error}"));
    assert_eq!(alignment.z_score, None);

    assert!(matches!(
        ce_align(
            &points,
            &points,
            CeOptions {
                fragment_similarity_threshold: -5.0,
                path_similarity_threshold: -4.0,
                ..CeOptions::original()
            }
        ),
        Err(CeError::InvalidOptions)
    ));
}

#[test]
fn biopython_ce_fixture_when_configured() {
    let (Some(reference), Some(mobile)) = (
        std::env::var_os("PDBIOX_CE_REFERENCE"),
        std::env::var_os("PDBIOX_CE_MOBILE"),
    ) else {
        return;
    };
    let reference = guide_coordinates(reference);
    let mobile = guide_coordinates(mobile);
    let alignment = ce_align(&reference, &mobile, CeOptions::original())
        .unwrap_or_else(|error| panic!("fixture CE failed: {error}"));
    assert!(
        (alignment.rmsd - 3.83).abs() < 0.02,
        "rmsd={}",
        alignment.rmsd
    );
}

fn guide_coordinates(path: std::ffi::OsString) -> Vec<[f32; 3]> {
    use pdbiox_core::io::{InputBuffer, ReadOptions};
    let bytes = std::fs::read(path).unwrap_or_else(|error| panic!("fixture read failed: {error}"));
    let input = InputBuffer::from_bytes(bytes);
    let (structure, _) = pdbiox_cif::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("fixture CIF failed: {findings:?}"));
    let mut chains: Vec<_> = structure.data().chains().collect();
    chains.sort_by(|left, right| left.label().cmp(&right.label()));
    let mut coordinates = Vec::new();
    for chain in chains {
        let mut residues: Vec<_> = chain.residues().collect();
        residues.sort_by_key(|residue| {
            residue
                .auth_seq_id()
                .or_else(|| residue.label_seq_id())
                .unwrap_or(i32::MIN)
        });
        coordinates.extend(residues.into_iter().filter_map(|residue| {
            residue
                .atom("CA")
                .or_else(|| residue.atom("C4'"))
                .and_then(pdbiox_core::structure::AtomRef::position)
        }));
    }
    coordinates
}
