//! In-process registration and catalog parity checks.

use pyo3::prelude::*;
use pyo3::types::PyModule;

#[cfg(feature = "full")]
type Scorer = fn(&[[f32; 3]], &[[f32; 3]]) -> Result<f64, molframe::compare::CompareError>;

#[test]
fn curated_module_and_catalog_are_bidirectionally_consistent() {
    Python::initialize();
    Python::attach(|py| {
        let module = match PyModule::new(py, "_native") {
            Ok(module) => module,
            Err(error) => panic!("test module should be created: {error}"),
        };
        if let Err(error) = crate::module::native(&module) {
            panic!("the extension module should register: {error}");
        }
        for capability in super::catalog::capabilities() {
            let domain = match module.getattr(capability.domain) {
                Ok(value) => value,
                Err(error) => panic!("{} should be registered: {error}", capability.domain),
            };
            if capability.eager {
                assert!(
                    domain.getattr(capability.name).is_ok(),
                    "{}.{} is catalogued but absent",
                    capability.domain,
                    capability.name
                );
            }
        }
        for removed in ["Plan", "Batch", "AtomIndex", "StructureData", "SpatialPlan"] {
            assert!(
                module.getattr(removed).is_err(),
                "{removed} leaked into root"
            );
        }
    });
}

#[cfg(feature = "full")]
#[test]
#[ignore = "requires NumPy on the embedded Python search path; run in the Python job"]
fn python_numeric_results_equal_rust_on_a_real_structure() {
    use numpy::{PyArray2, PyArrayMethods};
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../molframe-bench/data/1ubq.cif");
    let structure = molframe::read(path).expect("real fixture parses");
    let reference = structure.coordinates();
    let mobile: Vec<[f32; 3]> = reference
        .iter()
        .enumerate()
        .map(|(index, point)| {
            [
                point[0] + [0.0, 0.2, 0.4][index % 3],
                point[1] - 0.1,
                point[2] + 0.3,
            ]
        })
        .collect();
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "_native").expect("test module");
        crate::module::native(&module).expect("registers");
        let array = |points: &[[f32; 3]]| {
            PyArray2::from_vec2(
                py,
                &points.iter().map(|row| row.to_vec()).collect::<Vec<_>>(),
            )
            .expect("rectangular array")
        };
        let left = array(&mobile);
        let right = array(reference);
        let geometry = module.getattr("geometry").expect("geometry");
        let centroid: Option<[f64; 3]> = geometry
            .getattr("centroid")
            .expect("centroid")
            .call1((&right,))
            .expect("runs")
            .extract()
            .expect("point");
        assert_eq!(centroid, molframe::geometry::centroid(reference));
        let rmsd: f64 = geometry
            .getattr("rmsd")
            .expect("rmsd")
            .call1((&left, &right))
            .expect("runs")
            .extract()
            .expect("scalar");
        assert_eq!(
            rmsd.to_bits(),
            molframe::geometry::rmsd(&mobile, reference)
                .expect("native RMSD")
                .to_bits()
        );
        let compare = module.getattr("compare").expect("compare");
        let scorers: [(&str, Scorer); 3] = [
            ("tm_score", molframe::compare::tm_score),
            ("gdt_ts", molframe::compare::gdt_ts),
            ("gdt_ha", molframe::compare::gdt_ha),
        ];
        for (name, scorer) in scorers {
            let result: f64 = compare
                .getattr(name)
                .expect("registered scorer")
                .call1((&left, &right))
                .expect("runs")
                .extract()
                .expect("scalar");
            assert_eq!(
                result.to_bits(),
                scorer(&mobile, reference).expect("native score").to_bits(),
                "{name}"
            );
        }
        let matrix = geometry
            .getattr("distance_matrix")
            .expect("matrix")
            .call1((&right,))
            .expect("runs")
            .cast_into::<PyArray2<f64>>()
            .expect("matrix array");
        let native = molframe::geometry::distance_matrix(reference).expect("native matrix");
        assert_eq!(
            matrix.readonly().as_slice().expect("contiguous"),
            native.into_values()
        );
    });
}

#[cfg(feature = "full")]
#[test]
fn python_internal_coordinates_equal_the_native_forest_and_rebuild() {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../molframe-bench/data/1ubq.cif");
    let structure = molframe::read(path).expect("fixture");
    let native =
        molframe::ic::internal_coordinates(structure.engine(), molframe::ModelIndex::new(0))
            .expect("native forest");
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "_native").expect("module");
        crate::module::native(&module).expect("registers");
        let wrapped =
            Bound::new(py, crate::bindings::PyStructure::new(structure.clone())).expect("wrapper");
        let forest = module
            .getattr("geometry")
            .expect("geometry")
            .getattr("internal_coordinates")
            .expect("operation")
            .call1((&wrapped,))
            .expect("forest");
        let seeds: Vec<(usize, [f32; 3])> = forest
            .getattr("seeds")
            .expect("seeds")
            .extract()
            .expect("records");
        assert_eq!(
            seeds,
            native
                .seeds()
                .iter()
                .map(|(atom, point)| (atom.as_usize(), *point))
                .collect::<Vec<_>>()
        );
        let rebuilt: Vec<Option<[f32; 3]>> = forest
            .call_method0("rebuild")
            .expect("rebuild")
            .extract()
            .expect("positions");
        assert_eq!(rebuilt, native.rebuild().expect("native rebuild"));
        let positions: Vec<Option<[f32; 3]>> = structure
            .coordinates()
            .iter()
            .map(|point| Some([point[0] + 3.0, point[1] - 2.0, point[2] + 1.0]))
            .collect();
        let bat = forest
            .call_method1("measure_bat", (positions.clone(),))
            .expect("BAT");
        let measured: Vec<[f64; 3]> = bat
            .getattr("coordinates")
            .expect("coordinates")
            .extract()
            .expect("BAT rows");
        let expected = native.measure_bat(&positions).expect("native BAT");
        assert_eq!(measured, expected.coordinates());
        let rebuilt: Vec<Option<[f32; 3]>> = forest
            .call_method1("rebuild_bat", (&bat,))
            .expect("rebuild BAT")
            .extract()
            .expect("positions");
        assert_eq!(
            rebuilt,
            native.rebuild_bat(&expected).expect("native rebuilt BAT")
        );
    });
}

#[cfg(feature = "full")]
#[test]
#[ignore = "requires MOLFRAME_PARITY_CORPUS with the external Ligare refined-set"]
fn python_and_rust_read_the_same_external_ligand_corpus() {
    use numpy::{PyArray2, PyArrayMethods};
    let root = std::path::PathBuf::from(
        std::env::var("MOLFRAME_PARITY_CORPUS")
            .expect("set MOLFRAME_PARITY_CORPUS to the acquired refined-set"),
    );
    let mut paths = Vec::new();
    for directory in std::fs::read_dir(root).expect("external corpus") {
        let directory = directory.expect("directory").path();
        if directory.is_dir() {
            for entry in std::fs::read_dir(directory).expect("complex directory") {
                let path = entry.expect("entry").path();
                if path.extension().is_some_and(|extension| extension == "sdf") {
                    paths.push(path);
                }
            }
        }
    }
    paths.sort();
    assert!(
        !paths.is_empty(),
        "the external corpus must contain SDF inputs"
    );
    Python::initialize();
    Python::attach(|py| {
        crate::error::tests::install_package_for_tests(py);
        let module = PyModule::new(py, "_native").expect("module");
        crate::module::native(&module).expect("registers");
        let read = module.getattr("read").expect("read");
        for path in &paths {
            let native =
                molframe::read(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
            let wrapped = read
                .call1((path.to_string_lossy().as_ref(),))
                .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
            let coordinates = wrapped
                .getattr("coordinates")
                .expect("coordinates")
                .cast_into::<PyArray2<f32>>()
                .expect("coordinate array");
            let flat: Vec<f32> = native.coordinates().iter().flatten().copied().collect();
            assert_eq!(
                coordinates.readonly().as_slice().expect("contiguous"),
                flat,
                "{}",
                path.display()
            );
            let atoms: u32 = wrapped
                .getattr("atom_count")
                .expect("atom_count")
                .extract()
                .expect("count");
            assert_eq!(atoms, native.atom_count(), "{}", path.display());
        }
    });
    eprintln!(
        "Rust/Python coordinate and atom-count parity: {} external SDF inputs",
        paths.len()
    );
}
