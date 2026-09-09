//! Conversion of Rust diagnostics into Python exceptions.

use pdbiox::{Class, Code, DatasetError as NativeDatasetError, Diagnostic};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyTuple, PyType};

create_exception!(_native, PdbioxError, PyException);
create_exception!(_native, ParseError, PdbioxError);
create_exception!(_native, SchemaError, PdbioxError);
create_exception!(_native, ConsistencyError, PdbioxError);
create_exception!(_native, ConversionError, PdbioxError);
create_exception!(_native, GeometryError, PdbioxError);
create_exception!(_native, PolicyError, PdbioxError);
create_exception!(_native, PolicyConfigError, PdbioxError);
create_exception!(_native, DifferenceError, PdbioxError);
create_exception!(_native, ReexecutionError, PdbioxError);
create_exception!(_native, CapacityError, PdbioxError);
create_exception!(_native, TableError, PdbioxError);
create_exception!(_native, DatasetError, PdbioxError);
create_exception!(_native, ProviderError, PdbioxError);
create_exception!(_native, GraphError, PdbioxError);
create_exception!(_native, TableFileError, PdbioxError);
create_exception!(_native, DlpackError, PdbioxError);
create_exception!(_native, LoadError, PdbioxError);
create_exception!(_native, CompareError, PdbioxError);
create_exception!(_native, CadError, PdbioxError);
create_exception!(_native, CadConstructionError, PdbioxError);
create_exception!(_native, CeError, PdbioxError);
create_exception!(_native, GovernedCompareError, PdbioxError);
create_exception!(_native, BasePairError, PdbioxError);
create_exception!(_native, CationPiError, PdbioxError);
create_exception!(_native, DensityError, PdbioxError);
create_exception!(_native, DsspBinaryError, PdbioxError);
create_exception!(_native, DsspError, PdbioxError);
create_exception!(_native, DynamicsError, PdbioxError);
create_exception!(_native, EnsembleGeometryError, PdbioxError);
create_exception!(_native, EnsembleSimilarityError, PdbioxError);
create_exception!(_native, EnsembleStatisticsError, PdbioxError);
create_exception!(_native, FragmentMappingError, PdbioxError);
create_exception!(_native, GnmError, PdbioxError);
create_exception!(_native, GovernedAnalysisError, PdbioxError);
create_exception!(_native, GovernedEnsembleError, PdbioxError);
create_exception!(_native, GovernedNativeError, PdbioxError);
create_exception!(_native, HelicalError, PdbioxError);
create_exception!(_native, HseError, PdbioxError);
create_exception!(_native, HydrogenBondError, PdbioxError);
create_exception!(_native, KMeansError, PdbioxError);
create_exception!(_native, NativeError, PdbioxError);
create_exception!(_native, NucleicTorsionError, PdbioxError);
create_exception!(_native, PhysicalKernelError, PdbioxError);
create_exception!(_native, PiStackingError, PdbioxError);
create_exception!(_native, PolymerError, PdbioxError);
create_exception!(_native, PoreError, PdbioxError);
create_exception!(_native, RadialError, PdbioxError);
create_exception!(_native, StandaloneAnalysisError, PdbioxError);
create_exception!(_native, AtomDepthError, PdbioxError);
create_exception!(_native, BuriedSurfaceError, PdbioxError);
create_exception!(_native, SasaError, PdbioxError);
create_exception!(_native, SasaStreamError, PdbioxError);
create_exception!(_native, SurfaceGeometryError, PdbioxError);
create_exception!(_native, SurfaceWorkflowError, PdbioxError);
create_exception!(_native, MapStatisticsError, PdbioxError);
create_exception!(_native, MrcError, PdbioxError);
create_exception!(_native, MrcBrickError, PdbioxError);
create_exception!(_native, MonomerLibraryReadError, PdbioxError);
create_exception!(_native, RestraintError, PdbioxError);
create_exception!(_native, ReflectionError, PdbioxError);
create_exception!(_native, AltlocOccupancyError, PdbioxError);
create_exception!(_native, AltlocOccupancyKernelError, PdbioxError);
create_exception!(_native, BFactorError, PdbioxError);
create_exception!(_native, BFactorKernelError, PdbioxError);
create_exception!(_native, CcdCompletenessKernelError, PdbioxError);
create_exception!(_native, CompletenessError, PdbioxError);
create_exception!(_native, GovernedMapError, PdbioxError);
create_exception!(_native, LigandGeometryKernelError, PdbioxError);
create_exception!(_native, NucleicGeometryError, PdbioxError);
create_exception!(_native, PlanarityError, PdbioxError);
create_exception!(_native, PlaneRestraintKernelError, PdbioxError);
create_exception!(_native, RamachandranError, PdbioxError);
create_exception!(_native, RealSpaceCorrelationError, PdbioxError);
create_exception!(_native, ReferenceError, PdbioxError);
create_exception!(_native, RotamerError, PdbioxError);
create_exception!(_native, DielectricError, PdbioxError);
create_exception!(_native, DmsError, PdbioxError);
create_exception!(_native, ImdError, PdbioxError);
create_exception!(_native, MsdError, PdbioxError);
create_exception!(_native, PathSimilarityError, PdbioxError);
create_exception!(_native, TrajectoryError, PdbioxError);
create_exception!(_native, TrajectoryIoError, PdbioxError);
create_exception!(_native, UpdatingSelectionError, PdbioxError);
create_exception!(_native, WaterDynamicsError, PdbioxError);

static INDEX_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();
static KEY_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("PdbioxError", module.py().get_type::<PdbioxError>())?;
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
        mixed_exception(module.py(), "PdbioxIndexError", "IndexError")
    })?;
    let key_error = KEY_ERROR.get_or_try_init(module.py(), || {
        mixed_exception(module.py(), "PdbioxKeyError", "KeyError")
    })?;
    module.add("PdbioxIndexError", index_error.bind(module.py()))?;
    module.add("PdbioxKeyError", key_error.bind(module.py()))?;
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

pub(crate) fn graph_error(error: &pdbiox::GraphError) -> PyErr {
    GraphError::new_err(error.to_string())
}

pub(crate) fn table_file_error(error: &pdbiox::TableFileError) -> PyErr {
    TableFileError::new_err(error.to_string())
}

pub(crate) fn dlpack_error(error: &pdbiox::DlpackError) -> PyErr {
    DlpackError::new_err(error.to_string())
}

pub(crate) fn compare_error(error: &pdbiox::compare::CompareError) -> PyErr {
    CompareError::new_err(error.to_string())
}

pub(crate) fn cad_error(error: &pdbiox::compare::CadError) -> PyErr {
    CadError::new_err(error.to_string())
}

pub(crate) fn cad_construction_error(error: &pdbiox::compare::CadConstructionError) -> PyErr {
    CadConstructionError::new_err(error.to_string())
}

pub(crate) fn ce_error(error: &pdbiox::compare::CeError) -> PyErr {
    CeError::new_err(error.to_string())
}

pub(crate) fn governed_compare_error(error: &pdbiox::compare::GovernedCompareError) -> PyErr {
    GovernedCompareError::new_err(error.to_string())
}

fn mixed_exception(py: Python<'_>, name: &str, standard: &str) -> PyResult<Py<PyType>> {
    let builtins = py.import("builtins")?;
    let standard = builtins.getattr(standard)?.cast_into::<PyType>()?;
    let bases = PyTuple::new(
        py,
        [py.get_type::<PdbioxError>().into_any(), standard.into_any()],
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
        return PdbioxError::new_err(Code::E9001.cause());
    };
    match kind.bind(py).call1((value,)) {
        Ok(instance) => PyErr::from_value(instance),
        Err(error) => error,
    }
}

pub(crate) fn read_error(py: Python<'_>, findings: &[Diagnostic]) -> PyErr {
    let internal_invariant = Diagnostic::new(Code::E9001);
    let finding = match findings.first() {
        Some(finding) => finding,
        None => &internal_invariant,
    };
    let error = match finding.code().class() {
        Class::Syntax => ParseError::new_err(finding.message().to_owned()),
        Class::Schema => SchemaError::new_err(finding.message().to_owned()),
        Class::Consistency => ConsistencyError::new_err(finding.message().to_owned()),
        Class::Conversion => ConversionError::new_err(finding.message().to_owned()),
        Class::Geometry => GeometryError::new_err(finding.message().to_owned()),
        Class::Policy => PolicyError::new_err(finding.message().to_owned()),
        _ => PdbioxError::new_err(finding.message().to_owned()),
    };
    attach_diagnostic(py, error, finding)
}

pub(crate) fn cif_write_error(error: &pdbiox::CifWriteToError) -> PyErr {
    ConversionError::new_err(error.to_string())
}

fn attach_diagnostic(py: Python<'_>, error: PyErr, finding: &Diagnostic) -> PyErr {
    let value = error.value(py);
    let attributes = [
        ("code", finding.code().to_string()),
        ("message", finding.message().to_owned()),
        ("remedy", finding.remedy().to_owned()),
    ];
    for (name, attribute) in attributes {
        if let Err(set_error) = value.setattr(name, attribute) {
            return set_error;
        }
    }
    if let Some(span) = finding.span()
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
    error
}
