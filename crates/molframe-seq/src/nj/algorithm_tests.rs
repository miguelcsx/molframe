use super::neighbor_joining;

#[test]
fn additive_distances_recover_the_generating_tree() {
    // Distances from a tree with edges A-x=1, B-x=2, x-y=1, y-C=3, y-D=4.
    let labels = ["A", "B", "C", "D"];
    let distances = vec![
        vec![0.0, 3.0, 5.0, 6.0],
        vec![3.0, 0.0, 6.0, 7.0],
        vec![5.0, 6.0, 0.0, 7.0],
        vec![6.0, 7.0, 7.0, 0.0],
    ];
    let Some(tree) = neighbor_joining(&labels, &distances) else {
        panic!("valid matrix");
    };
    // A and B join first with their true branch lengths of 1 and 2.
    assert_eq!(tree.to_newick(), "(((A:1,B:2):1,C:3):0,D:4);");
}

#[test]
fn two_taxa_join_directly() {
    let Some(tree) = neighbor_joining(&["A", "B"], &[vec![0.0, 4.0], vec![4.0, 0.0]]) else {
        panic!("valid");
    };
    assert_eq!(tree.to_newick(), "(A:0,B:4);");
}

#[test]
fn a_non_square_matrix_is_rejected() {
    assert!(neighbor_joining(&["A", "B"], &[vec![0.0]]).is_none());
}

mod validation {
    use super::super::{
        NegativeBranches, NjError, NjOptions, neighbor_joining_with, validate_distance_matrix,
    };

    fn three(matrix: [[f64; 3]; 3]) -> Vec<Vec<f64>> {
        matrix.iter().map(|row| row.to_vec()).collect()
    }

    const LABELS: [&str; 3] = ["A", "B", "C"];

    #[test]
    fn rejects_non_finite_negative_asymmetric_and_nonzero_diagonal() {
        let nan = three([[0.0, f64::NAN, 1.0], [f64::NAN, 0.0, 1.0], [1.0, 1.0, 0.0]]);
        assert!(matches!(
            validate_distance_matrix(&LABELS, &nan),
            Err(NjError::NotFinite { .. })
        ));
        let negative = three([[0.0, -1.0, 1.0], [-1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]);
        assert!(matches!(
            validate_distance_matrix(&LABELS, &negative),
            Err(NjError::Negative { .. })
        ));
        let asymmetric = three([[0.0, 1.0, 1.0], [2.0, 0.0, 1.0], [1.0, 1.0, 0.0]]);
        assert!(matches!(
            validate_distance_matrix(&LABELS, &asymmetric),
            Err(NjError::Asymmetric { row: 0, column: 1 })
        ));
        let diagonal = three([[0.5, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]);
        assert!(matches!(
            validate_distance_matrix(&LABELS, &diagonal),
            Err(NjError::NonZeroDiagonal { index: 0 })
        ));
        assert_eq!(validate_distance_matrix(&[], &[]), Err(NjError::Empty));
    }

    #[test]
    fn negative_branches_follow_the_policy() {
        // A non-additive but valid matrix that yields a negative branch.
        let labels = ["A", "B", "C", "D"];
        let distances = vec![
            vec![0.0, 1.0, 10.0, 10.0],
            vec![1.0, 0.0, 1.0, 10.0],
            vec![10.0, 1.0, 0.0, 1.0],
            vec![10.0, 10.0, 1.0, 0.0],
        ];
        let keep = neighbor_joining_with(&labels, &distances, NjOptions::default());
        let clamp = neighbor_joining_with(
            &labels,
            &distances,
            NjOptions {
                negative_branches: NegativeBranches::ClampToZero,
            },
        );
        let (Ok(keep), Ok(clamp)) = (keep, clamp) else {
            panic!("valid matrix");
        };
        assert!(keep.to_newick().contains(":-"), "{}", keep.to_newick());
        assert!(!clamp.to_newick().contains(":-"), "{}", clamp.to_newick());
    }
}
