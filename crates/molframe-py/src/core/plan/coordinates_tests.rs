use super::*;
use numpy::IntoPyArray;
use numpy::ndarray::arr2;

#[test]
fn coordinate_comparisons_validate_pair_shape_and_alignment_at_construction() {
    Python::initialize();
    Python::attach(|py| {
        if py.import("numpy").is_err() {
            return;
        }
        assert_message(
            PyRmsd::new(py, coordinates(py, 2), coordinates(py, 1)),
            "mobile and reference must contain the same number of coordinates",
        );
        let invalid_shape = arr2(&[[0.0_f32, 0.0]]).into_pyarray(py).unbind().into_any();
        assert_message(
            PyLddt::new(py, invalid_shape, coordinates(py, 1), 15.0),
            "coordinates must have shape (n, 3)",
        );
        if let Err(error) = PyTmScore::new(py, coordinates(py, 2), coordinates(py, 2)) {
            panic!("aligned coordinate buffers should construct: {error}");
        }
    });
}

fn coordinates(py: Python<'_>, count: usize) -> Py<PyAny> {
    let values = match count {
        1 => arr2(&[[0.0_f32, 0.0, 0.0]]),
        2 => arr2(&[[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        _ => panic!("test only defines one or two coordinates"),
    };
    values.into_pyarray(py).unbind().into_any()
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
