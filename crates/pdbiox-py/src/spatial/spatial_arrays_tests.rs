use super::*;
use numpy::ndarray::{Array1, Array2, arr1, arr2};
use numpy::{
    IntoPyArray, PyArray1, PyArray2, PyArrayDescrMethods, PyArrayMethods, PyReadonlyArray1,
    PyReadonlyArray2, PyUntypedArrayMethods, dtype,
};
use pyo3::types::{PyAnyMethods, PyModule, PySlice};

fn positions(py: Python<'_>) -> PyReadonlyArray2<'_, f32> {
    arr2(&[
        [0.0_f32, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
    ])
    .into_pyarray(py)
    .readonly()
}

fn indices(py: Python<'_>) -> PyReadonlyArray1<'_, u32> {
    arr1(&[0_u32, 1, 2, 3]).into_pyarray(py).readonly()
}

fn numpy_available(py: Python<'_>) -> bool {
    py.import("numpy").is_ok()
}

fn table_indices<'py>(table: &Bound<'py, PyNeighborTable>) -> Bound<'py, PyArray2<u32>> {
    match table
        .getattr("indices")
        .and_then(|array| Ok(array.cast_into::<PyArray2<u32>>()?))
    {
        Ok(array) => array,
        Err(error) => panic!("neighbor table must expose an index array: {error}"),
    }
}

fn table_distances<'py>(table: &Bound<'py, PyNeighborTable>) -> Bound<'py, PyArray1<f32>> {
    match table
        .getattr("distance")
        .and_then(|array| Ok(array.cast_into::<PyArray1<f32>>()?))
    {
        Ok(array) => array,
        Err(error) => panic!("neighbor table must expose a distance array: {error}"),
    }
}

#[test]
fn array_selections_must_already_be_sorted_and_unique() {
    let Err(error) = array_selection(Some(&[2, 0, 2]), 3) else {
        panic!("unsorted array selection must be rejected");
    };
    assert_eq!(error, "atom index arrays must be sorted and unique");
}

#[test]
fn neighbor_table_columns_are_zero_copy_read_only_and_owner_retained() {
    Python::initialize();
    Python::attach(|py| {
        if !numpy_available(py) {
            return;
        }
        let result = match neighbor_pairs(
            py,
            positions(py),
            1.1,
            Some(indices(py)),
            None,
            PySpatialBackend::BruteForce,
            None,
        ) {
            Ok(result) => result,
            Err(error) => panic!("neighbor query should succeed: {error}"),
        };
        let table = match Py::new(py, result) {
            Ok(table) => table,
            Err(error) => panic!("neighbor table should become a Python object: {error}"),
        };
        let expected_indices = table.bind(py).borrow().indices.as_ptr();
        let expected_distances = table.bind(py).borrow().distance.as_ptr();
        let pairs = table_indices(table.bind(py));
        let distances = table_distances(table.bind(py));

        assert!(pairs.dtype().is_equiv_to(&dtype::<u32>(py)));
        assert!(distances.dtype().is_equiv_to(&dtype::<f32>(py)));
        assert_eq!(pairs.shape(), [1, 2]);
        assert_eq!(distances.shape(), [1]);
        assert_eq!(pairs.data().cast_const(), expected_indices);
        assert_eq!(distances.data().cast_const(), expected_distances);
        assert!(pairs.try_readwrite().is_err());
        assert!(distances.try_readwrite().is_err());

        drop(table);
        let pair_view = pairs.readonly();
        let pair_values = match pair_view.as_slice() {
            Ok(values) => values,
            Err(error) => panic!("neighbor index view should stay contiguous: {error}"),
        };
        let distance_view = distances.readonly();
        let distance_values = match distance_view.as_slice() {
            Ok(values) => values,
            Err(error) => panic!("neighbor distance view should stay contiguous: {error}"),
        };
        assert_eq!(pair_values, [0, 1]);
        assert_eq!(distance_values, [1.0]);
    });
}

#[test]
fn canonical_queries_return_columnar_or_index_arrays() {
    Python::initialize();
    Python::attach(|py| {
        if !numpy_available(py) {
            return;
        }
        let profile = PySpatialSearchOptions::balanced_value();
        let pairs = match neighbor_pairs_with_options(
            py,
            positions(py),
            1.1,
            profile,
            Some(indices(py)),
            None,
            None,
        ) {
            Ok(result) => result,
            Err(error) => panic!("profiled neighbor query should succeed: {error}"),
        };
        assert_eq!(pairs.distance, [1.0]);
        assert_eq!(pairs.indices, [0, 1]);

        let within = match atoms_within(
            py,
            positions(py),
            arr1(&[0_u32]).into_pyarray(py).readonly(),
            1.1,
            Some(indices(py)),
            PySpatialBackend::BruteForce,
            None,
        ) {
            Ok(result) => result,
            Err(error) => panic!("within query should succeed: {error}"),
        };
        assert!(within.dtype().is_equiv_to(&dtype::<u32>(py)));
        assert_eq!(within.shape(), [2]);
        assert!(within.try_readwrite().is_err());
        let within_view = within.readonly();
        let within_values = match within_view.as_slice() {
            Ok(values) => values,
            Err(error) => panic!("within result should be contiguous: {error}"),
        };
        assert_eq!(within_values, [0, 1]);

        let profiled_within = match atoms_within_with_options(
            py,
            positions(py),
            arr1(&[0_u32]).into_pyarray(py).readonly(),
            1.1,
            profile,
            Some(indices(py)),
            None,
        ) {
            Ok(result) => result,
            Err(error) => panic!("profiled within query should succeed: {error}"),
        };
        let profiled_view = profiled_within.readonly();
        let profiled_values = match profiled_view.as_slice() {
            Ok(values) => values,
            Err(error) => panic!("profiled within result should be contiguous: {error}"),
        };
        assert_eq!(profiled_values, [0, 1]);

        let nearest = match nearest_neighbors(py, positions(py), indices(py), 0, 2, None, None) {
            Ok(result) => result,
            Err(error) => panic!("nearest-neighbor query should succeed: {error}"),
        };
        assert_eq!(nearest.indices, [0, 1, 0, 2]);
        assert_eq!(nearest.distance, [1.0, 3.0]);
    });
}

#[test]
fn array_bindings_reject_wrong_dtypes_and_coordinate_shapes() {
    Python::initialize();
    Python::attach(|py| {
        if !numpy_available(py) {
            return;
        }
        let floating = Array2::<f64>::zeros((2, 3)).into_pyarray(py);
        let integer = Array1::<i64>::zeros(2).into_pyarray(py);
        assert!(floating.extract::<PyReadonlyArray2<'_, f32>>().is_err());
        assert!(integer.extract::<PyReadonlyArray1<'_, u32>>().is_err());

        let invalid_positions = arr2(&[[0.0_f32, 0.0], [1.0, 0.0]])
            .into_pyarray(py)
            .readonly();
        assert!(
            neighbor_pairs(
                py,
                invalid_positions,
                1.1,
                None,
                None,
                PySpatialBackend::BruteForce,
                None,
            )
            .is_err()
        );
    });
}

#[test]
fn canonical_bindings_reject_noncontiguous_coordinates_and_indices() {
    Python::initialize();
    Python::attach(|py| {
        if !numpy_available(py) {
            return;
        }
        let coordinates =
            arr2(&[[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]).into_pyarray(py);
        let noncontiguous_coordinates = match coordinates.transpose() {
            Ok(array) => array.readonly(),
            Err(error) => panic!("transpose should create a view: {error}"),
        };
        assert!(
            neighbor_pairs(
                py,
                noncontiguous_coordinates,
                1.1,
                None,
                None,
                PySpatialBackend::BruteForce,
                None,
            )
            .is_err()
        );

        let source = arr1(&[0_u32, 99, 1, 99, 2, 99]).into_pyarray(py);
        let view = match source.get_item(PySlice::new(py, 0, 6, 2)) {
            Ok(view) => view,
            Err(error) => panic!("stride slice should create a view: {error}"),
        };
        let noncontiguous_indices = match view.extract::<PyReadonlyArray1<'_, u32>>() {
            Ok(indices) => indices,
            Err(error) => panic!("stride slice should remain a uint32 array: {error}"),
        };
        assert!(
            atoms_within(
                py,
                positions(py),
                noncontiguous_indices,
                1.1,
                None,
                PySpatialBackend::BruteForce,
                None,
            )
            .is_err()
        );
    });
}

#[test]
fn canonical_functions_are_registered_without_array_or_alias_paths() {
    Python::initialize();
    Python::attach(|py| {
        let module = match PyModule::new(py, "spatial") {
            Ok(module) => module,
            Err(error) => panic!("test module should be created: {error}"),
        };
        if let Err(error) = super::super::register(&module) {
            panic!("spatial functions should register: {error}");
        }
        for name in [
            "NeighborTable",
            "neighbor_pairs",
            "neighbor_pairs_with_options",
            "atoms_within",
            "atoms_within_with_options",
            "nearest_neighbors",
        ] {
            assert!(module.getattr(name).is_ok(), "missing {name}");
        }
        for name in [
            "neighbor_pairs_array",
            "neighbor_pairs_array_with_options",
            "atoms_within_array",
            "atoms_within_array_with_options",
            "nearest_neighbors_array",
            "pairs_within",
            "pairs_within_with_options",
            "within",
            "within_with_options",
        ] {
            assert!(
                module.getattr(name).is_err(),
                "obsolete {name} is registered"
            );
        }
    });
}
