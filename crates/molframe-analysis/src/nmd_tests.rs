use super::{NmdError, read_nmd, write_nmd};

const NMD: &str = "name test modes
atomnames CA CA CA
resnames ALA GLY SER
chainids A A A
resids 1 2 3
bfactors 10.0 11.0 12.0
coordinates 0.0 0.0 0.0 3.8 0.0 0.0 7.6 0.0 0.0
mode 1 12.5 0.1 0.0 0.0 0.0 0.2 0.0 0.0 0.0 0.3
mode 2 8.0 0.0 0.1 0.0 0.1 0.0 0.0 0.0 0.0 0.1
";

#[test]
fn a_mode_file_reads_its_sites_labels_and_modes() {
    let Ok(set) = read_nmd(NMD.as_bytes()) else {
        panic!("the fixture file reads")
    };
    assert_eq!(set.name.as_deref(), Some("test modes"));
    assert_eq!(set.len(), 3);
    assert_eq!(set.residue_ids, vec![1, 2, 3]);
    assert_eq!(set.residue_names.len(), 3);
    assert_eq!(set.chain_ids.len(), 3);
    assert_eq!(set.modes.len(), 2);
    assert_eq!(set.modes[0].index, 1);
    assert!((set.modes[0].scale - 12.5).abs() < 1e-9);
    assert_eq!(set.modes[0].displacements.len(), 3);
    assert!((set.modes[0].displacements[2][2] - 0.3).abs() < 1e-9);
}

#[test]
fn displacing_applies_every_amplitude_to_the_reference_coordinates() {
    let Ok(set) = read_nmd(NMD.as_bytes()) else {
        panic!("the fixture file reads")
    };
    let mut moved = vec![[0.0_f64; 3]; set.len()];
    set.displace(&[2.0, 0.0], &mut moved);
    assert!((moved[0][0] - 0.2).abs() < 1e-9, "first site was {moved:?}");
    assert!((moved[1][1] - 0.4).abs() < 1e-9);
    assert!(
        (moved[2][0] - 7.6).abs() < 1e-9,
        "unmoved axes are preserved"
    );
}

#[test]
fn a_mode_covering_the_wrong_number_of_sites_is_refused() {
    let short = NMD.replace(
        "mode 1 12.5 0.1 0.0 0.0 0.0 0.2 0.0 0.0 0.0 0.3",
        "mode 1 12.5 0.1 0.0 0.0",
    );
    assert!(matches!(
        read_nmd(short.as_bytes()),
        Err(NmdError::SiteMismatch { expected: 3, .. })
    ));
}

#[test]
fn a_mode_before_the_coordinate_block_is_refused() {
    let reordered = "mode 1 1.0 0.1 0.0 0.0\ncoordinates 0.0 0.0 0.0\n";
    assert!(matches!(
        read_nmd(reordered.as_bytes()),
        Err(NmdError::ModeBeforeCoordinates)
    ));
}

#[test]
fn a_coordinate_block_that_is_not_whole_vectors_is_refused() {
    let ragged = NMD.replace(
        "coordinates 0.0 0.0 0.0 3.8 0.0 0.0 7.6 0.0 0.0",
        "coordinates 0.0 0.0 0.0 3.8",
    );
    assert!(matches!(
        read_nmd(ragged.as_bytes()),
        Err(NmdError::UnalignedVectors {
            keyword: "coordinates"
        })
    ));
}

#[test]
fn a_mode_written_without_its_index_and_scale_still_reads() {
    let bare = "coordinates 0.0 0.0 0.0 1.0 0.0 0.0\nmode 0.1 0.0 0.0 0.2 0.0 0.0\n";
    let Ok(set) = read_nmd(bare.as_bytes()) else {
        panic!("a bare mode record reads")
    };
    assert_eq!(set.modes.len(), 1);
    assert_eq!(set.modes[0].index, 0);
    assert!((set.modes[0].scale - 1.0).abs() < 1e-9);
    assert!((set.modes[0].displacements[1][0] - 0.2).abs() < 1e-9);
}

#[test]
fn writing_and_reading_a_mode_set_returns_the_same_model() {
    let Ok(set) = read_nmd(NMD.as_bytes()) else {
        panic!("the fixture file reads")
    };
    let written = write_nmd(&set);
    let Ok(reread) = read_nmd(written.as_bytes()) else {
        panic!("the written file reads back")
    };
    assert_eq!(set, reread);
}

#[test]
fn unknown_viewer_records_are_ignored_rather_than_refused() {
    let styled = format!("nmwiz_load something.nmd\n{NMD}scale 3.0\n");
    let Ok(set) = read_nmd(styled.as_bytes()) else {
        panic!("styling records do not stop the read")
    };
    assert_eq!(set.modes.len(), 2);
}
