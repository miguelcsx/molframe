use super::*;
use molframe_core::execution::TempStoragePolicy;

fn context(root: &Path) -> ExecutionContext {
    ExecutionContext::builder()
        .temp_storage_policy(TempStoragePolicy::directory(root, 10_000_000))
        .build()
        .expect("context")
}

#[test]
fn invalid_inputs_are_refused_before_starting_an_executable_or_reserving_storage() {
    let context = ExecutionContext::default();
    let atoms = SurfaceAtoms {
        positions: &[[0.0; 3]],
        radii: &[1.7],
        rows: &[1],
    };
    assert!(matches!(
        native_molecular_surface(
            Path::new("missing"),
            NativeSurfaceProvider::Surf,
            atoms,
            NativeSurfaceOptions::default(),
            &context
        ),
        Err(NativeSurfaceError::Invalid(_))
    ));
    assert_eq!(context.spill_bytes(), 0);
}

#[cfg(unix)]
#[test]
fn unsuccessful_providers_retain_bounded_diagnostics_and_remove_private_files() {
    use std::os::unix::fs::PermissionsExt;
    let root = tempfile::tempdir().expect("root");
    let executable = root.path().join("provider");
    std::fs::write(
        &executable,
        "#!/bin/sh\ni=0\nwhile [ $i -lt 1500 ]; do echo intermediate-diagnostic; i=$((i+1)); done\necho 'ERROR Too many RS components' >&2\nexit 9\n",
    )
    .expect("script");
    std::fs::set_permissions(&executable, std::fs::Permissions::from_mode(0o700))
        .expect("permissions");
    let context = context(root.path());
    let error = native_molecular_surface(
        &executable,
        NativeSurfaceProvider::Msms,
        SurfaceAtoms {
            positions: &[[0.0; 3]],
            radii: &[1.7],
            rows: &[0],
        },
        NativeSurfaceOptions {
            temporary_bytes: 10_000_000,
            ..NativeSurfaceOptions::default()
        },
        &context,
    )
    .expect_err("provider fails");
    let NativeSurfaceError::Failed { status, diagnostic } = error else {
        panic!("native status and diagnostics survive the process failure");
    };
    assert_eq!(status.code(), Some(9));
    assert!(diagnostic.len() <= 8192);
    assert!(diagnostic.ends_with("ERROR Too many RS components"));
    assert_eq!(context.spill_bytes(), 0);
    assert_eq!(context.reserved_bytes(), 0);
    assert_eq!(std::fs::read_dir(root.path()).expect("root").count(), 1);
}

#[cfg(unix)]
#[test]
fn successful_status_without_complete_triangulation_is_a_provider_output_error() {
    use std::os::unix::fs::PermissionsExt;
    let root = tempfile::tempdir().expect("root");
    let executable = root.path().join("provider");
    let context = context(root.path());
    for script in [
        "#!/bin/sh\nexit 0\n",
        "#!/bin/sh\ntouch surface.vert\nexit 0\n",
    ] {
        std::fs::write(&executable, script).expect("provider");
        std::fs::set_permissions(&executable, std::fs::Permissions::from_mode(0o700))
            .expect("permissions");
        assert!(matches!(
            native_molecular_surface(
                &executable,
                NativeSurfaceProvider::Msms,
                SurfaceAtoms {
                    positions: &[[0.0; 3]],
                    radii: &[1.7],
                    rows: &[0]
                },
                NativeSurfaceOptions {
                    temporary_bytes: 10_000_000,
                    ..NativeSurfaceOptions::default()
                },
                &context,
            ),
            Err(NativeSurfaceError::Output(SurfaceOutputError::Invalid(_)))
        ));
        assert_eq!(context.reserved_bytes(), 0);
        assert_eq!(context.spill_bytes(), 0);
        assert_eq!(std::fs::read_dir(root.path()).expect("root").count(), 1);
    }
}

#[test]
#[ignore = "requires caller-installed scientific providers"]
fn actual_providers_return_complete_analytic_sphere_geometry_with_original_atom_identity() {
    let root = tempfile::tempdir().expect("root");
    let context = context(root.path());
    let atoms = SurfaceAtoms {
        positions: &[[99.0; 3], [0.0; 3]],
        radii: &[1.7; 2],
        rows: &[1],
    };
    for (provider, variable) in [
        (NativeSurfaceProvider::Surf, "MOLFRAME_SURF"),
        (NativeSurfaceProvider::Msms, "MOLFRAME_MSMS"),
    ] {
        let executable = std::env::var_os(variable).expect("provider executable");
        let surface = native_molecular_surface(
            Path::new(&executable),
            provider,
            atoms,
            NativeSurfaceOptions {
                temporary_bytes: 10_000_000,
                ..NativeSurfaceOptions::default()
            },
            &context,
        )
        .expect("scientific geometry");
        assert!(surface.mesh.geometry.faces.len() >= 64);
        assert!(surface.mesh.atom_rows.iter().all(|&row| row == 1));
        for position in &surface.mesh.geometry.vertices {
            let radius = position
                .iter()
                .map(|value| value * value)
                .sum::<f32>()
                .sqrt();
            assert!((radius - 1.7).abs() < 0.02, "{provider:?}: {radius}");
        }
        assert!(context.reserved_bytes() > 0);
        assert_eq!(context.spill_bytes(), 0);
        drop(surface);
        assert_eq!(context.reserved_bytes(), 0);
        assert_eq!(
            std::fs::read_dir(root.path()).expect("directory").count(),
            0
        );
    }
}

#[cfg(unix)]
#[test]
fn cancellation_and_timeout_reap_the_child_before_removing_its_private_files() {
    use std::os::unix::fs::PermissionsExt;
    let root = tempfile::tempdir().expect("root");
    let executable = root.path().join("provider");
    std::fs::write(&executable, "#!/bin/sh\nexec /bin/sleep 10\n").expect("script");
    std::fs::set_permissions(&executable, std::fs::Permissions::from_mode(0o700))
        .expect("permissions");
    let context = context(root.path());
    let atoms = SurfaceAtoms {
        positions: &[[0.0; 3]],
        radii: &[1.7],
        rows: &[0],
    };
    let options = NativeSurfaceOptions {
        timeout: Duration::from_millis(40),
        temporary_bytes: 10_000_000,
        ..NativeSurfaceOptions::default()
    };
    assert!(matches!(
        native_molecular_surface(
            &executable,
            NativeSurfaceProvider::Surf,
            atoms,
            options,
            &context
        ),
        Err(NativeSurfaceError::Timeout)
    ));
    assert_eq!(context.spill_bytes(), 0);
    let token = context.cancellation().clone();
    std::thread::scope(|scope| {
        scope.spawn(move || {
            std::thread::sleep(Duration::from_millis(30));
            token.cancel();
        });
        assert!(matches!(
            native_molecular_surface(
                &executable,
                NativeSurfaceProvider::Surf,
                atoms,
                options,
                &context
            ),
            Err(NativeSurfaceError::Cancelled(_))
        ));
    });
    assert_eq!(context.spill_bytes(), 0);
    assert_eq!(
        std::fs::read_dir(root.path()).expect("directory").count(),
        1
    );
}

#[test]
#[ignore = "requires caller-installed MSMS"]
fn all_disconnected_msms_components_are_retained_with_their_original_atom_rows() {
    let root = tempfile::tempdir().expect("root");
    let context = context(root.path());
    let executable = std::env::var_os("MOLFRAME_MSMS").expect("provider executable");
    let atoms = SurfaceAtoms {
        positions: &[[0.0; 3], [10.0, 0.0, 0.0]],
        radii: &[1.7; 2],
        rows: &[0, 1],
    };
    let options = NativeSurfaceOptions {
        temporary_bytes: 10_000_000,
        ..NativeSurfaceOptions::default()
    };
    let surface = native_molecular_surface(
        Path::new(&executable),
        NativeSurfaceProvider::Msms,
        atoms,
        options,
        &context,
    )
    .expect("all components");
    assert!(surface.mesh.atom_rows.contains(&0));
    assert!(surface.mesh.atom_rows.contains(&1));
    assert_eq!(surface.mesh.geometry.report.components, 2);
    assert_eq!(surface.mesh.geometry.report.boundary_edges, 0);
}

#[test]
fn isolated_component_provenance_checks_obey_task_cancellation() {
    let context = ExecutionContext::default();
    context.cancellation().cancel();
    let mut mesh = MolecularSurfaceMesh {
        geometry: crate::IndexedSurfaceMesh {
            vertices: vec![[1.7, 0.0, 0.0]],
            vertex_normals: vec![[1.0, 0.0, 0.0]],
            ..crate::IndexedSurfaceMesh::default()
        },
        atom_rows: vec![0].into_boxed_slice(),
        discarded_degenerate_faces: 0,
    };
    let atoms = SurfaceAtoms {
        positions: &[[0.0; 3]],
        radii: &[1.7],
        rows: &[0],
    };
    assert!(matches!(
        provenance::validate_isolated_sphere(
            &mut mesh,
            NativeAtoms {
                atoms,
                sample: None
            },
            &context,
        ),
        Err(NativeSurfaceError::Cancelled(_))
    ));
}
