use super::*;
use numpy::IntoPyArray;
use numpy::ndarray::{arr1, arr2};

#[test]
fn spatial_operations_validate_cutoffs_at_construction() {
    Python::initialize();
    Python::attach(|py| {
        if py.import("numpy").is_err() {
            return;
        }
        assert_message(
            PyNeighborPairs::new(py, positions(py), f32::NAN, None, None, None, None),
            "cutoff must be finite and non-negative",
        );
        assert_message(
            PyAtomsWithin::new(py, positions(py), targets(py), -0.1, None, None, None),
            "cutoff must be finite and non-negative",
        );
        if let Err(error) = PyNeighborPairs::new(py, positions(py), 0.0, None, None, None, None) {
            panic!("zero cutoff is a valid native spatial search: {error}");
        }
    });
}

fn positions(py: Python<'_>) -> Py<PyAny> {
    arr2(&[[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0]])
        .into_pyarray(py)
        .unbind()
        .into_any()
}

fn targets(py: Python<'_>) -> Py<PyAny> {
    arr1(&[0_u32]).into_pyarray(py).unbind().into_any()
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
