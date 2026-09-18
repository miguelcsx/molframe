use super::*;
use numpy::IntoPyArray;
use numpy::ndarray::{Array1, arr2};

#[test]
fn surface_operations_require_atom_aligned_radii_and_masks() {
    Python::initialize();
    Python::attach(|py| {
        if py.import("numpy").is_err() {
            return;
        }
        let short_radii = Array1::from_vec(vec![1.5_f32])
            .into_pyarray(py)
            .unbind()
            .into_any();
        assert_message(
            PySasa::new(py, positions(py), short_radii, 1.4, 96),
            "radii must contain one value per coordinate",
        );

        let radii = Array1::from_vec(vec![1.5_f32, 1.7])
            .into_pyarray(py)
            .unbind()
            .into_any();
        let short_mask = Array1::from_vec(vec![true])
            .into_pyarray(py)
            .unbind()
            .into_any();
        assert_message(
            PyBuriedSurfaceOperation::new(py, positions(py), radii, short_mask, 1.4, 96),
            "first must contain one mask value per coordinate",
        );

        let valid_radii = Array1::from_vec(vec![1.5_f32, 1.7])
            .into_pyarray(py)
            .unbind()
            .into_any();
        if let Err(error) = PySasa::new(py, positions(py), valid_radii, 1.4, 96) {
            panic!("atom-aligned surface inputs should construct: {error}");
        }
    });
}

fn positions(py: Python<'_>) -> Py<PyAny> {
    arr2(&[[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0]])
        .into_pyarray(py)
        .unbind()
        .into_any()
}

fn assert_message<T>(result: PyResult<T>, expected: &str) {
    let Err(error) = result else {
        panic!("operation construction should fail");
    };
    assert!(
        error.to_string().contains(expected),
        "expected {expected:?} in {error}"
    );
}
