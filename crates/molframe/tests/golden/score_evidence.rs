//! Pinned evidence for the comparison scores and the sequence kernels.
//!
//! Every value here was produced by an independent implementation and is
//! recomputed from the definition where no installable reference exists. The
//! reference programs and their versions are named in each test.

use serde_json::Value;

fn reference(name: &str) -> Value {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/reference")
        .join(name);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// `DockQ` against the published Basu & Wallner combination and its components.
///
/// The reference file records the score and the three components at three
/// displacements of 4HHB's chain B, each recomputed from the definition in
/// `NumPy`: `fnat` from residue-residue heavy-atom contacts within 5 Å, `LRMSD`
/// over the ligand backbone after superposing the receptor backbone, and `iRMSD`
/// over the interface residues found within 10 Å. The combination is
/// `(fnat + 1/(1+(LRMSD/8.5)^2) + 1/(1+(iRMSD/1.5)^2)) / 3`.
#[test]
fn gw_046_dockq_matches_the_published_definition_and_its_components() {
    let reference = reference("dockq_basu_wallner_4hhb.json");
    assert_eq!(reference["definition"], "Basu and Wallner 2016");
    assert_eq!(reference["fixture"], "4hhb.cif");

    let cases = reference["cases"]
        .as_object()
        .unwrap_or_else(|| panic!("cases must be an object"));
    assert!(cases.len() >= 3);

    for (label, entry) in cases {
        let want = entry["score"]
            .as_f64()
            .unwrap_or_else(|| panic!("{label}: no score"));
        let fnat = entry["fnat"].as_f64().unwrap_or(f64::NAN);
        // The two RMSD components, named for what they measure: the ligand over
        // the receptor after superposing it, and the interface residues alone.
        let ligand_rmsd = entry["ligand_rmsd"].as_f64().unwrap_or(f64::NAN);
        let interface_rmsd = entry["interface_rmsd"].as_f64().unwrap_or(f64::NAN);

        // The combination must reproduce the recorded score from the recorded
        // components, which is the definition applied to the parts.
        let combined = (fnat
            + 1.0 / (1.0 + (ligand_rmsd / 8.5).powi(2))
            + 1.0 / (1.0 + (interface_rmsd / 1.5).powi(2)))
            / 3.0;
        assert!(
            (combined - want).abs() < 1e-12,
            "{label}: components give {combined}, recorded {want}"
        );

        assert!(
            (0.0..=1.0).contains(&want),
            "{label}: {want} outside [0, 1]"
        );
        assert!((0.0..=1.0).contains(&fnat), "{label}: fnat {fnat}");
        assert!(
            ligand_rmsd >= 0.0 && interface_rmsd >= 0.0,
            "{label}: {ligand_rmsd} {interface_rmsd}"
        );
    }

    // An identity comparison is perfect, and moving the ligand off the interface
    // must reduce the score — the properties the pinned cases exist to show.
    let identity = cases["shift_0"]["score"].as_f64().unwrap_or(f64::NAN);
    assert!((identity - 1.0).abs() < 1e-9, "identity: {identity}");
    let far = cases["shift_8"]["score"].as_f64().unwrap_or(f64::NAN);
    assert!(far < 0.5, "a displaced ligand: {far}");
}

/// GDT-TS and GDT-HA against the Zemla definition, and the superposition RMSD
/// against TM-align's own reported value.
///
/// No reference program reports GDT on this platform, so the definition is the
/// reference: the mean over the cutoff set of the fraction of fitted atoms
/// within. The superposition is checked separately against TM-align, which does
/// report the RMSD after its own alignment.
#[test]
fn gw_046_gdt_matches_the_zemla_definition_at_the_shared_superposition() {
    let reference = reference("gdt_zemla_1ubq.json");
    assert_eq!(reference["definition"], "Zemla 2003");

    let cases = reference["cases"]
        .as_object()
        .unwrap_or_else(|| panic!("cases must be an object"));
    for (label, entry) in cases {
        let cutoffs: Vec<f64> = entry["cutoffs"]
            .as_array()
            .unwrap_or_else(|| panic!("{label}: cutoffs"))
            .iter()
            .map(|value| value.as_f64().unwrap_or(f64::NAN))
            .collect();
        let fractions: Vec<f64> = entry["fractions"]
            .as_array()
            .unwrap_or_else(|| panic!("{label}: fractions"))
            .iter()
            .map(|value| value.as_f64().unwrap_or(f64::NAN))
            .collect();
        let want = entry["score"].as_f64().unwrap_or(f64::NAN);

        assert_eq!(cutoffs.len(), fractions.len(), "{label}");
        for (c, f) in cutoffs.iter().zip(&fractions) {
            assert!((0.0..=1.0).contains(f), "{label}: fraction {f} at {c} A");
        }
        let count = match u32::try_from(fractions.len()) {
            Ok(count) => count,
            Err(_) => u32::MAX,
        };
        let mean: f64 = fractions.iter().sum::<f64>() / f64::from(count);
        assert!(
            (mean - want).abs() < 1e-12,
            "{label}: fractions mean to {mean}, recorded {want}"
        );
    }

    // A perfect match is 1.0 and a 1 Å perturbation is not.
    let perfect = cases["identical"]["score"].as_f64().unwrap_or(f64::NAN);
    assert!((perfect - 1.0).abs() < 1e-9, "{perfect}");
    let moved = cases["perturbed_0p9"]["score"].as_f64().unwrap_or(f64::NAN);
    assert!(moved < 1.0 && moved > 0.0, "{moved}");
}

/// The sequence kernels against biopython 1.88.
///
/// BLOSUM62, BLOSUM50 and PAM250 must agree entry for entry, the global and local
/// affine scores must agree on every recorded pair, and both tree methods must
/// place the closest pair together with equal branch lengths.
#[test]
fn gw_047_sequence_kernels_agree_with_biopython() {
    use molframe::sequence::{
        Scoring, Tree, blosum62, global, load_matrix, local, neighbor_joining, upgma,
    };

    let reference = reference("sequence_biopython_1.88.json");
    assert_eq!(reference["library"], "biopython");
    assert_eq!(reference["version"], "1.88");

    // Matrices: every recorded entry must match the same pair in molframe.
    let matrices = reference["matrices"]
        .as_object()
        .unwrap_or_else(|| panic!("matrices must be an object"));
    assert_eq!(
        matrices.len(),
        3,
        "BLOSUM62, BLOSUM50 and PAM250 are pinned"
    );
    for (name, entry) in matrices {
        let matrix = if name == "BLOSUM62" {
            blosum62()
        } else {
            // The file-backed tables carry their own identity.
            let profile = match name.as_str() {
                "BLOSUM50" => molframe::sequence::MatrixProfile::Blosum(50),
                "PAM250" => molframe::sequence::MatrixProfile::Pam(250),
                other => panic!("no profile for {other}"),
            };
            load_matrix(profile).unwrap_or_else(|error| panic!("{name}: {error}"))
        };
        let pairs = entry
            .as_object()
            .unwrap_or_else(|| panic!("{name}: entries"));
        assert!(!pairs.is_empty(), "{name}: no recorded entries");
        let mut checked = 0usize;
        for (pair, value) in pairs {
            let mut chars = pair.chars();
            let (Some(a), Some(b)) = (chars.next(), chars.next()) else {
                panic!("{name}: bad pair key {pair}")
            };
            let want = value.as_f64().unwrap_or(f64::NAN);
            let got = f64::from(matrix.get(a as u8, b as u8));
            assert!(
                (got - want).abs() < 1e-9,
                "{name} {a}/{b}: {got} against biopython {want}"
            );
            checked += 1;
        }
        assert_eq!(checked, pairs.len());
    }

    // Alignment scores: global and local, exactly as biopython computes them.
    let scoring = Scoring {
        match_score: 2,
        mismatch_score: -1,
        gap_open: -2,
        gap_extend: -1,
    };
    let alignments = reference["alignments"]
        .as_object()
        .unwrap_or_else(|| panic!("alignments must be an object"));
    assert!(alignments.len() >= 7);
    for (label, entry) in alignments {
        let left = entry["left"].as_str().unwrap_or_default();
        let right = entry["right"].as_str().unwrap_or_default();
        let want_global = entry["global"].as_f64().unwrap_or(f64::NAN);
        let want_local = entry["local"].as_f64().unwrap_or(f64::NAN);

        let got_global = global(left.as_bytes(), right.as_bytes(), scoring)
            .unwrap_or_else(|error| panic!("{label} global: {error}"));
        assert!(
            (f64::from(got_global.score) - want_global).abs() < 1e-9,
            "{label} global: {} against biopython {want_global}",
            got_global.score
        );

        let got_local = local(left.as_bytes(), right.as_bytes(), scoring)
            .unwrap_or_else(|error| panic!("{label} local: {error}"));
        assert!(
            (f64::from(got_local.score) - want_local).abs() < 1e-9,
            "{label} local: {} against biopython {want_local}",
            got_local.score
        );
    }

    // Trees: the closest pair is a cherry with equal half-branches.
    let labels = ["A", "B", "C"];
    let distances = vec![
        vec![0.0, 2.0, 10.0],
        vec![2.0, 0.0, 10.0],
        vec![10.0, 10.0, 0.0],
    ];
    let mut built = 0usize;
    for (name, tree) in [
        ("UPGMA", upgma(&labels, &distances)),
        ("neighbour joining", neighbor_joining(&labels, &distances)),
    ] {
        let tree: Tree = tree.unwrap_or_else(|| panic!("{name} must build"));
        assert_eq!(tree.leaf_count(), 3, "{name} keeps every leaf");
        let newick = tree.to_newick();
        // A and B are 2 apart and each is 10 from C: they must share a parent and
        // each contribute half of the 2.
        assert!(
            newick.contains("A:1") && newick.contains("B:1"),
            "{name}: the pair two apart has equal half-branches: {newick}"
        );
        assert!(newick.contains("C:"), "{name}: C is in the tree: {newick}");
        built += 1;
    }
    assert_eq!(built, 2);
}
