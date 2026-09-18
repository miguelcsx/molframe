use super::*;
use numpy::IntoPyArray;
use numpy::ndarray::Array3;

#[test]
fn trajectory_constructors_reject_invalid_controls_before_native_execution() {
    Python::initialize();
    Python::attach(|py| {
        if py.import("numpy").is_err() {
            return;
        }
        assert_message(
            PyRmsdToReference::new(py, frames(py), 2, Some(PyFrameAlignment::Unaligned), None),
            "reference must index one of 2 trajectory frames",
        );
        assert_message(
            PyMeanSquaredDisplacementOperation::new(py, frames(py), 2, None, None),
            "maximum_lag must be smaller than the frame count",
        );
        assert_message(
            PyPairwiseFittedRmsd::new(py, frames(py), Some(0), None),
            "memory_limit must be positive",
        );
        assert_message(
            PyGeneralizedProcrustesMean::new(py, frames(py), Some(f64::NAN), Some(10), None),
            "tolerance must be finite and positive",
        );
        assert_message(
            PyGeneralizedProcrustesMean::new(py, frames(py), Some(1.0e-6), Some(0), None),
            "maximum_iterations must be positive",
        );
        if let Err(error) =
            PyRmsdToReference::new(py, frames(py), 1, Some(PyFrameAlignment::Unaligned), None)
        {
            panic!("valid trajectory controls should construct: {error}");
        }
    });
}

fn frames(py: Python<'_>) -> Py<PyAny> {
    Array3::<f32>::zeros((2, 2, 3))
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
