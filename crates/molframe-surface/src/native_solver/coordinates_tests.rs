use super::super::{
    NativeSurfaceOptions, NativeSurfaceProvider, native_molecular_surface_at_sample,
};
use super::*;
use molframe_core::ExecutionContext;
use std::path::Path;

#[test]
fn surface_input_samples_only_selected_rows_without_materializing_a_frame() {
    let start = [[-100.0; 3], [2.0, 4.0, 6.0]];
    let end = [[100.0; 3], [6.0, 8.0, 10.0]];
    let atoms = NativeAtoms {
        atoms: SurfaceAtoms {
            positions: &start,
            radii: &[1.0, 1.7],
            rows: &[1],
        },
        sample: Some(SurfaceCoordinateSample {
            end_positions: &end,
            weight: 0.25,
        }),
    };
    let context = ExecutionContext::builder()
        .temp_storage_policy(molframe_core::TempStoragePolicy::directory(
            std::env::temp_dir(),
            10000,
        ))
        .build()
        .expect("bounded files");
    super::super::validate(atoms, NativeSurfaceOptions::default()).expect("sample");
    let workspace = context
        .temporary_workspace_with_file_limit(10000, 2)
        .expect("files");
    for (provider, expected) in [
        (NativeSurfaceProvider::Surf, "0 1.7 3 5 7\n"),
        (NativeSurfaceProvider::Msms, "3 5 7 1.7\n"),
    ] {
        super::super::process::write_input(workspace.path(), provider, atoms, &context)
            .expect("input");
        assert_eq!(
            std::fs::read_to_string(workspace.path().join("atoms")).expect("file"),
            expected
        );
    }
    assert_eq!(atoms.atoms.positions.as_ptr(), start.as_ptr());
    assert_eq!(
        atoms.sample.expect("sample").end_positions.as_ptr(),
        end.as_ptr()
    );
}

#[test]
fn invalid_samples_fail_before_launching_a_provider_or_reserving_temporary_files() {
    let context = ExecutionContext::default();
    let atoms = SurfaceAtoms {
        positions: &[[0.0; 3]],
        radii: &[1.7],
        rows: &[0],
    };
    for sample in [
        SurfaceCoordinateSample {
            end_positions: &[],
            weight: 0.5,
        },
        SurfaceCoordinateSample {
            end_positions: &[[0.0; 3]],
            weight: -0.1,
        },
        SurfaceCoordinateSample {
            end_positions: &[[0.0; 3]],
            weight: 1.1,
        },
        SurfaceCoordinateSample {
            end_positions: &[[0.0; 3]],
            weight: f32::NAN,
        },
        SurfaceCoordinateSample {
            end_positions: &[[f32::INFINITY; 3]],
            weight: 0.0,
        },
    ] {
        assert!(matches!(
            native_molecular_surface_at_sample(
                Path::new("absent"),
                NativeSurfaceProvider::Surf,
                atoms,
                sample,
                NativeSurfaceOptions::default(),
                &context
            ),
            Err(NativeSurfaceError::Invalid(_))
        ));
    }
    assert_eq!(context.spill_bytes(), 0);
}
