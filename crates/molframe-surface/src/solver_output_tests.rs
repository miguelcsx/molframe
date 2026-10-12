use super::*;

const VERTICES: &str = "0 0 0 0 0 2 1 1 1\n1 0 0 0 0 2 1 2 1\n0 1 0 0 0 2 1 3 1\n";
const FACES: &str = "1 2 3 1 1\n";
const SURF: &str = "0\n0 0 0 0 0 2\n1 0 0 0 0 2\n0 1 0 0 0 2\n";

#[test]
fn surf_shares_only_identical_positions_analytic_normals_and_original_atom_owners() {
    let neighbor = "0\n0 0 0 0 0 2\n0 1 0 0 0 2\n-1 0 0 0 0 2\n";
    let text = format!("{SURF}{neighbor}");
    let cancellation = CancellationToken::new();
    let surface = read_surf_surface(
        text.as_bytes(),
        &[73],
        SurfaceOutputLimits {
            vertices: 4,
            ..SurfaceOutputLimits::default()
        },
        &cancellation,
    )
    .expect("shared analytical vertices");
    assert_eq!(surface.geometry.vertices.len(), 4);
    assert_eq!(
        surface.geometry.faces,
        [SurfaceFace([0, 1, 2]), SurfaceFace([0, 2, 3])]
    );
    assert_eq!(surface.geometry.report.components, 1);
    assert_eq!(surface.geometry.report.boundary_edges, 4);
    assert_eq!(surface.geometry.area().to_bits(), 1.0f64.to_bits());
    for changed in [
        neighbor.replacen("0\n", "1\n", 1),
        neighbor.replace("0 0 2", "2 0 0"),
    ] {
        let text = format!("{SURF}{changed}");
        let surface = read_surf_surface(
            text.as_bytes(),
            &[73, 91],
            SurfaceOutputLimits::default(),
            &cancellation,
        )
        .expect("attribute seams");
        assert_eq!(surface.geometry.vertices.len(), 6);
        assert_eq!(surface.geometry.area().to_bits(), 1.0f64.to_bits());
    }
}

#[test]
fn collapsed_surf_faces_are_counted_and_removed_without_changing_valid_triangles() {
    let collapsed = "0\n0 0 0 0 0 2\n1 0 0 0 0 2\n2 0 0 0 0 2\n";
    let text = format!("{collapsed}{SURF}");
    let surface = read_surf_surface(
        text.as_bytes(),
        &[73],
        SurfaceOutputLimits::default(),
        &CancellationToken::new(),
    )
    .expect("native tessellation");
    assert_eq!(surface.discarded_degenerate_faces, 1);
    assert_eq!(surface.geometry.vertices.len(), 3);
    assert_eq!(surface.geometry.faces, [SurfaceFace([0, 1, 2])]);
    assert_eq!(surface.geometry.area().to_bits(), 0.5f64.to_bits());
    assert!(
        read_surf_surface(
            collapsed.as_bytes(),
            &[73],
            SurfaceOutputLimits::default(),
            &CancellationToken::new()
        )
        .is_err()
    );
    let limits = SurfaceOutputLimits {
        faces: 1,
        ..SurfaceOutputLimits::default()
    };
    assert!(matches!(
        read_surf_surface(text.as_bytes(), &[73], limits, &CancellationToken::new()),
        Err(SurfaceOutputError::Limit("faces"))
    ));
}

#[test]
fn msms_preserves_source_atom_rows_and_normalizes_analytic_normals() {
    let surface = read_msms_surface(
        VERTICES.as_bytes(),
        FACES.as_bytes(),
        &[9, 17, 25],
        SurfaceOutputLimits::default(),
        &CancellationToken::new(),
    )
    .expect("MSMS surface");
    assert_eq!(&*surface.atom_rows, &[9, 17, 25]);
    assert_eq!(surface.geometry.faces, [SurfaceFace([0, 1, 2])]);
    for normal in &surface.geometry.vertex_normals {
        assert_eq!(normal.map(f32::to_bits), [0.0, 0.0, 1.0].map(f32::to_bits));
    }
    assert_eq!(surface.geometry.area().to_bits(), 0.5f64.to_bits());
}

#[test]
fn collapsed_msms_faces_are_counted_without_discarding_the_provider_vertex_columns() {
    let vertices = format!("{VERTICES}2 0 0 0 0 2 1 1 1\n");
    let faces = format!("1 2 4 1 1\n{FACES}");
    let cancellation = CancellationToken::new();
    let surface = read_msms_surface(
        vertices.as_bytes(),
        faces.as_bytes(),
        &[9, 17, 25],
        SurfaceOutputLimits::default(),
        &cancellation,
    )
    .expect("rounded native triangles");
    assert_eq!(surface.discarded_degenerate_faces, 1);
    assert_eq!(surface.geometry.vertices.len(), 4);
    assert_eq!(&*surface.atom_rows, &[9, 17, 25, 9]);
    assert_eq!(surface.geometry.faces, [SurfaceFace([0, 1, 2])]);
    assert_eq!(surface.geometry.area().to_bits(), 0.5f64.to_bits());
    assert!(matches!(
        read_msms_surface(
            vertices.as_bytes(),
            faces.as_bytes(),
            &[9, 17, 25],
            SurfaceOutputLimits {
                faces: 1,
                ..SurfaceOutputLimits::default()
            },
            &cancellation
        ),
        Err(SurfaceOutputError::Limit("faces"))
    ));
}

#[test]
fn surf_records_may_span_lines_without_losing_their_zero_based_atom_map() {
    for text in [SURF.to_owned(), SURF.replace('\n', " ")] {
        let surface = read_surf_surface(
            text.as_bytes(),
            &[73],
            SurfaceOutputLimits::default(),
            &CancellationToken::new(),
        )
        .expect("Surf surface");
        assert_eq!(&*surface.atom_rows, &[73; 3]);
        assert_eq!(surface.geometry.faces, [SurfaceFace([0, 1, 2])]);
        assert_eq!(surface.geometry.area().to_bits(), 0.5f64.to_bits());
    }
}

#[test]
fn invalid_msms_records_never_return_an_incomplete_prefix() {
    for vertices in [
        VERTICES.replace("0 0 2", "0 0 0"),
        VERTICES.replace("1 3 1", "1 4 1"),
        VERTICES.replace("1 1 1", "1 0 1"),
        VERTICES.replace("1 0 0", "NaN 0 0"),
        format!("{VERTICES}0 0\n"),
        format!("{VERTICES}{}\n", "0 ".repeat(600)),
    ] {
        assert!(
            read_msms_surface(
                vertices.as_bytes(),
                FACES.as_bytes(),
                &[9, 17, 25],
                SurfaceOutputLimits::default(),
                &CancellationToken::new()
            )
            .is_err()
        );
    }
    for faces in [
        "1 2 4 1 1\n",
        "1 1 2 1 1\n",
        "0 2 3 1 1\n",
        "1 2 3 1\n",
        "1 2 3 1 1 7\n",
    ] {
        assert!(
            read_msms_surface(
                VERTICES.as_bytes(),
                faces.as_bytes(),
                &[9, 17, 25],
                SurfaceOutputLimits::default(),
                &CancellationToken::new()
            )
            .is_err()
        );
    }
}

#[test]
fn surface_readers_enforce_vertex_face_byte_and_cancellation_limits() {
    for limits in [
        SurfaceOutputLimits {
            vertices: 2,
            ..SurfaceOutputLimits::default()
        },
        SurfaceOutputLimits {
            faces: 0,
            ..SurfaceOutputLimits::default()
        },
        SurfaceOutputLimits {
            bytes: 5,
            ..SurfaceOutputLimits::default()
        },
    ] {
        assert!(matches!(
            read_msms_surface(
                VERTICES.as_bytes(),
                FACES.as_bytes(),
                &[9, 17, 25],
                limits,
                &CancellationToken::new()
            ),
            Err(SurfaceOutputError::Limit(_))
        ));
        assert!(matches!(
            read_surf_surface(SURF.as_bytes(), &[73], limits, &CancellationToken::new()),
            Err(SurfaceOutputError::Limit(_))
        ));
    }
    let cancellation = CancellationToken::new();
    cancellation.cancel();
    assert!(matches!(
        read_surf_surface(
            SURF.as_bytes(),
            &[73],
            SurfaceOutputLimits::default(),
            &cancellation
        ),
        Err(SurfaceOutputError::Cancelled(_))
    ));
    assert!(
        read_surf_surface(
            b"0\n0 0".as_slice(),
            &[73],
            SurfaceOutputLimits::default(),
            &CancellationToken::new()
        )
        .is_err()
    );
}
