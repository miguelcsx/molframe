//! Conversion of Rust diagnostics into Python exceptions.

use crate::contract::PyDiagnostic;
use molframe::{Class, Code, DatasetError as NativeDatasetError, Diagnostic};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyTuple, PyType};

create_exception!(_native, MolframeError, PyException);
create_exception!(_native, ParseError, MolframeError);
create_exception!(_native, SchemaError, MolframeError);
create_exception!(_native, ConsistencyError, MolframeError);
create_exception!(_native, ConversionError, MolframeError);
create_exception!(_native, GeometryError, MolframeError);
create_exception!(_native, PolicyError, MolframeError);
create_exception!(_native, PolicyConfigError, MolframeError);
create_exception!(_native, DifferenceError, MolframeError);
create_exception!(_native, ReexecutionError, MolframeError);
create_exception!(_native, CapacityError, MolframeError);
create_exception!(_native, TableError, MolframeError);
create_exception!(_native, DatasetError, MolframeError);
create_exception!(_native, ProviderError, MolframeError);
create_exception!(_native, GraphError, MolframeError);
create_exception!(_native, TableFileError, MolframeError);
create_exception!(_native, DlpackError, MolframeError);
create_exception!(_native, LoadError, MolframeError);
create_exception!(_native, CompareError, MolframeError);
create_exception!(_native, CadError, MolframeError);
create_exception!(_native, CadConstructionError, MolframeError);
create_exception!(_native, CeError, MolframeError);
create_exception!(_native, GovernedCompareError, MolframeError);
create_exception!(_native, BasePairError, MolframeError);
create_exception!(_native, CationPiError, MolframeError);
create_exception!(_native, DensityError, MolframeError);
create_exception!(_native, DsspBinaryError, MolframeError);
create_exception!(_native, DsspError, MolframeError);
create_exception!(_native, DynamicsError, MolframeError);
create_exception!(_native, EnsembleGeometryError, MolframeError);
create_exception!(_native, EnsembleSimilarityError, MolframeError);
create_exception!(_native, EnsembleStatisticsError, MolframeError);
create_exception!(_native, FragmentMappingError, MolframeError);
create_exception!(_native, GnmError, MolframeError);
create_exception!(_native, GovernedAnalysisError, MolframeError);
create_exception!(_native, GovernedEnsembleError, MolframeError);
create_exception!(_native, GovernedNativeError, MolframeError);
create_exception!(_native, HelicalError, MolframeError);
create_exception!(_native, HseError, MolframeError);
create_exception!(_native, HydrogenBondError, MolframeError);
create_exception!(_native, KMeansError, MolframeError);
create_exception!(_native, NativeError, MolframeError);
create_exception!(_native, NucleicTorsionError, MolframeError);
create_exception!(_native, PhysicalKernelError, MolframeError);
create_exception!(_native, PiStackingError, MolframeError);
create_exception!(_native, PolymerError, MolframeError);
create_exception!(_native, PoreError, MolframeError);
create_exception!(_native, RadialError, MolframeError);
create_exception!(_native, StandaloneAnalysisError, MolframeError);
create_exception!(_native, AtomDepthError, MolframeError);
create_exception!(_native, BuriedSurfaceError, MolframeError);
create_exception!(_native, SasaError, MolframeError);
create_exception!(_native, SasaStreamError, MolframeError);
create_exception!(_native, SurfaceGeometryError, MolframeError);
create_exception!(_native, SurfaceWorkflowError, MolframeError);
create_exception!(_native, MapStatisticsError, MolframeError);
create_exception!(_native, MrcError, MolframeError);
create_exception!(_native, MrcBrickError, MolframeError);
create_exception!(_native, MonomerLibraryReadError, MolframeError);
create_exception!(_native, RestraintError, MolframeError);
create_exception!(_native, ReflectionError, MolframeError);
create_exception!(_native, AltlocOccupancyError, MolframeError);
create_exception!(_native, AltlocOccupancyKernelError, MolframeError);
create_exception!(_native, BFactorError, MolframeError);
create_exception!(_native, BFactorKernelError, MolframeError);
create_exception!(_native, CcdCompletenessKernelError, MolframeError);
create_exception!(_native, CompletenessError, MolframeError);
create_exception!(_native, GovernedMapError, MolframeError);
create_exception!(_native, LigandGeometryKernelError, MolframeError);
create_exception!(_native, NucleicGeometryError, MolframeError);
create_exception!(_native, PlanarityError, MolframeError);
create_exception!(_native, PlaneRestraintKernelError, MolframeError);
create_exception!(_native, RamachandranError, MolframeError);
create_exception!(_native, RealSpaceCorrelationError, MolframeError);
create_exception!(_native, ReferenceError, MolframeError);
create_exception!(_native, RotamerError, MolframeError);
create_exception!(_native, DielectricError, MolframeError);
create_exception!(_native, DmsError, MolframeError);
create_exception!(_native, ImdError, MolframeError);
create_exception!(_native, MsdError, MolframeError);
create_exception!(_native, PathSimilarityError, MolframeError);
create_exception!(_native, TrajectoryError, MolframeError);
create_exception!(_native, TrajectoryIoError, MolframeError);
create_exception!(_native, UpdatingSelectionError, MolframeError);
create_exception!(_native, WaterDynamicsError, MolframeError);

static INDEX_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();
static KEY_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("MolframeError", module.py().get_type::<MolframeError>())?;
    module.add("ParseError", module.py().get_type::<ParseError>())?;
    module.add("SchemaError", module.py().get_type::<SchemaError>())?;
    module.add(
        "ConsistencyError",
        module.py().get_type::<ConsistencyError>(),
    )?;
    module.add("ConversionError", module.py().get_type::<ConversionError>())?;
    module.add("GeometryError", module.py().get_type::<GeometryError>())?;
    module.add("PolicyError", module.py().get_type::<PolicyError>())?;
    module.add(
        "PolicyConfigError",
        module.py().get_type::<PolicyConfigError>(),
    )?;
    module.add("DifferenceError", module.py().get_type::<DifferenceError>())?;
    module.add(
        "ReexecutionError",
        module.py().get_type::<ReexecutionError>(),
    )?;
    module.add("CapacityError", module.py().get_type::<CapacityError>())?;
    module.add("TableError", module.py().get_type::<TableError>())?;
    module.add("DatasetError", module.py().get_type::<DatasetError>())?;
    module.add("ProviderError", module.py().get_type::<ProviderError>())?;
    module.add("GraphError", module.py().get_type::<GraphError>())?;
    module.add("TableFileError", module.py().get_type::<TableFileError>())?;
    module.add("DlpackError", module.py().get_type::<DlpackError>())?;
    module.add("LoadError", module.py().get_type::<LoadError>())?;
    module.add("CompareError", module.py().get_type::<CompareError>())?;
    module.add("CadError", module.py().get_type::<CadError>())?;
    module.add(
        "CadConstructionError",
        module.py().get_type::<CadConstructionError>(),
    )?;
    module.add("CeError", module.py().get_type::<CeError>())?;
    module.add(
        "GovernedCompareError",
        module.py().get_type::<GovernedCompareError>(),
    )?;
    register_analysis_errors(module)?;
    register_trajectory_errors(module)?;
    let index_error = INDEX_ERROR.get_or_try_init(module.py(), || {
        mixed_exception(module.py(), "MolframeIndexError", "IndexError")
    })?;
    let key_error = KEY_ERROR.get_or_try_init(module.py(), || {
        mixed_exception(module.py(), "MolframeKeyError", "KeyError")
    })?;
    module.add("MolframeIndexError", index_error.bind(module.py()))?;
    module.add("MolframeKeyError", key_error.bind(module.py()))?;
    Ok(())
}

fn register_trajectory_errors(module: &Bound<'_, PyModule>) -> PyResult<()> {
    macro_rules! add {
        ($($name:ident),+ $(,)?) => {
            $(module.add(stringify!($name), module.py().get_type::<$name>())?;)+
        };
    }
    add!(
        DielectricError,
        DmsError,
        ImdError,
        MsdError,
        PathSimilarityError,
        TrajectoryError,
        TrajectoryIoError,
        UpdatingSelectionError,
        WaterDynamicsError,
    );
    Ok(())
}

fn register_analysis_errors(module: &Bound<'_, PyModule>) -> PyResult<()> {
    macro_rules! add {
        ($($name:ident),+ $(,)?) => {
            $(module.add(stringify!($name), module.py().get_type::<$name>())?;)+
        };
    }
    add!(
        BasePairError,
        CationPiError,
        DensityError,
        DsspBinaryError,
        DsspError,
        DynamicsError,
        EnsembleGeometryError,
        EnsembleSimilarityError,
        EnsembleStatisticsError,
        FragmentMappingError,
        GnmError,
        GovernedAnalysisError,
        GovernedEnsembleError,
        GovernedNativeError,
        HelicalError,
        HseError,
        HydrogenBondError,
        KMeansError,
        NativeError,
        NucleicTorsionError,
        PhysicalKernelError,
        PiStackingError,
        PolymerError,
        PoreError,
        RadialError,
        StandaloneAnalysisError,
        AtomDepthError,
        BuriedSurfaceError,
        SasaError,
        SasaStreamError,
        SurfaceGeometryError,
        SurfaceWorkflowError,
        MapStatisticsError,
        MrcError,
        MrcBrickError,
        MonomerLibraryReadError,
        RestraintError,
        ReflectionError,
        AltlocOccupancyError,
        AltlocOccupancyKernelError,
        BFactorError,
        BFactorKernelError,
        CcdCompletenessKernelError,
        CompletenessError,
        GovernedMapError,
        LigandGeometryKernelError,
        NucleicGeometryError,
        PlanarityError,
        PlaneRestraintKernelError,
        RamachandranError,
        RealSpaceCorrelationError,
        ReferenceError,
        RotamerError,
    );
    Ok(())
}

pub(crate) fn dataset_error(error: NativeDatasetError) -> PyErr {
    match error {
        NativeDatasetError::IndexOutOfBounds { index } => index_error_for_dataset(index),
        error => DatasetError::new_err(error.to_string()),
    }
}

fn index_error_for_dataset(index: usize) -> PyErr {
    pyo3::exceptions::PyIndexError::new_err(index)
}

pub(crate) fn graph_error(error: &molframe::GraphError) -> PyErr {
    GraphError::new_err(error.to_string())
}

pub(crate) fn table_file_error(error: &molframe::TableFileError) -> PyErr {
    TableFileError::new_err(error.to_string())
}

pub(crate) fn dlpack_error(error: &molframe::DlpackError) -> PyErr {
    DlpackError::new_err(error.to_string())
}

pub(crate) fn compare_error(error: &molframe::compare::CompareError) -> PyErr {
    CompareError::new_err(error.to_string())
}

pub(crate) fn cad_error(error: &molframe::compare::CadError) -> PyErr {
    CadError::new_err(error.to_string())
}

pub(crate) fn cad_construction_error(error: &molframe::compare::CadConstructionError) -> PyErr {
    CadConstructionError::new_err(error.to_string())
}

pub(crate) fn ce_error(error: &molframe::compare::CeError) -> PyErr {
    CeError::new_err(error.to_string())
}

pub(crate) fn governed_compare_error(error: &molframe::compare::GovernedCompareError) -> PyErr {
    GovernedCompareError::new_err(error.to_string())
}

fn mixed_exception(py: Python<'_>, name: &str, standard: &str) -> PyResult<Py<PyType>> {
    let builtins = py.import("builtins")?;
    let standard = builtins.getattr(standard)?.cast_into::<PyType>()?;
    let bases = PyTuple::new(
        py,
        [
            py.get_type::<MolframeError>().into_any(),
            standard.into_any(),
        ],
    )?;
    let namespace = PyDict::new(py);
    Ok(builtins
        .getattr("type")?
        .call1((name, bases, namespace))?
        .cast_into::<PyType>()
        .map(Bound::unbind)?)
}

pub(crate) fn index_error(py: Python<'_>, index: isize) -> PyErr {
    instantiate_standard(py, &INDEX_ERROR, index)
}

pub(crate) fn key_error(py: Python<'_>, key: &str) -> PyErr {
    instantiate_standard(py, &KEY_ERROR, key)
}

fn instantiate_standard(
    py: Python<'_>,
    kind: &PyOnceLock<Py<PyType>>,
    value: impl for<'a> IntoPyObject<'a>,
) -> PyErr {
    let Some(kind) = kind.get(py) else {
        return MolframeError::new_err(Code::E9001.cause());
    };
    match kind.bind(py).call1((value,)) {
        Ok(instance) => PyErr::from_value(instance),
        Err(error) => error,
    }
}

pub(crate) fn read_error(py: Python<'_>, findings: &[Diagnostic]) -> PyErr {
    let internal_invariant = Diagnostic::new(Code::E9001);
    let head = match findings.first() {
        Some(head) => head,
        None => &internal_invariant,
    };
    let error = match head.code().class() {
        Class::Syntax => ParseError::new_err(head.message().to_owned()),
        Class::Schema => SchemaError::new_err(head.message().to_owned()),
        Class::Consistency => ConsistencyError::new_err(head.message().to_owned()),
        Class::Conversion => ConversionError::new_err(head.message().to_owned()),
        Class::Geometry => GeometryError::new_err(head.message().to_owned()),
        Class::Policy => PolicyError::new_err(head.message().to_owned()),
        _ => MolframeError::new_err(head.message().to_owned()),
    };
    attach_diagnostic(py, error, findings)
}

/// Attaches the finding set to a raised error.
///
/// The scalar attributes describe the *head* of the list, which is the worst
/// finding because the core orders them worst-first — so a caller that wants the
/// reason can read a string. `findings` carries the whole ordered list, because a
/// reader reports every problem it noticed and an exception that named only the
/// first would lose the rest.
fn attach_diagnostic(py: Python<'_>, error: PyErr, findings: &[Diagnostic]) -> PyErr {
    let internal_invariant = Diagnostic::new(Code::E9001);
    let head = match findings.first() {
        Some(head) => head,
        None => &internal_invariant,
    };
    let value = error.value(py);
    let attributes = [
        ("code", head.code().to_string()),
        ("message", head.message().to_owned()),
        ("remedy", head.remedy().to_owned()),
    ];
    for (name, attribute) in attributes {
        if let Err(set_error) = value.setattr(name, attribute) {
            return set_error;
        }
    }
    if let Some(span) = head.span()
        && let Err(set_error) = value.setattr(
            "span",
            (
                span.start.byte_offset,
                span.end,
                span.start.line,
                span.start.column,
            ),
        )
    {
        return set_error;
    }
    let attached: Vec<PyDiagnostic> = findings.iter().map(PyDiagnostic::from).collect();
    if let Err(set_error) = value.setattr("findings", attached) {
        return set_error;
    }
    error
}
