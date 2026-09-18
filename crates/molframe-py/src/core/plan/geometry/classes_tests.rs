use super::*;
use numpy::IntoPyArray;
use numpy::ndarray::{Array1, Array3, arr2};

#[test]
fn geometry_constructors_validate_borrowed_array_contracts_before_execution() {
    Python::initialize();
    Python::attach(|py| {
        if py.import("numpy").is_err() {
            return;
        }
        let invalid_shape = arr2(&[[0.0_f32, 0.0], [1.0, 0.0]])
            .into_pyarray(py)
            .unbind()
            .into_any();
        assert_message(
            PyCentroid::new(py, invalid_shape),
            "coordinates must have shape (n, 3)",
        );

        let positions = coordinates(py);
        let masses = Array1::from_vec(vec![12.0_f64])
            .into_pyarray(py)
            .unbind()
            .into_any();
        assert_message(
            PyCentreOfMass::new(py, positions, Some(masses)),
            "masses must contain one value per coordinate",
        );

        let invalid_frames = Array3::<f64>::zeros((2, 2, 3))
            .into_pyarray(py)
            .unbind()
            .into_any();
        assert!(PyRmsf::new(py, invalid_frames).is_err());

        let valid_positions = coordinates(py);
        let valid_masses = Array1::from_vec(vec![12.0_f64, 16.0])
            .into_pyarray(py)
            .unbind()
            .into_any();
        if let Err(error) = PyCentreOfMass::new(py, valid_positions, Some(valid_masses)) {
            panic!("aligned contiguous geometry inputs should construct: {error}");
        }
    });
}

fn coordinates(py: Python<'_>) -> Py<PyAny> {
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
