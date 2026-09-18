//! The names most programs want in scope.
//!
//! `use molframe::prelude::*;` brings in the verbs ([`read`], [`write`]), the
//! vocabulary they speak (a [`Structure`], its [`ReadOptions`], the [`Diagnostic`]
//! that reports what was wrong with a file), and the pieces needed to call them.
//! Everything is gated by the feature that provides it, so a caller linking one
//! format gets exactly that format's names.

// The data, and what a read says about it.
pub use crate::{
    Analysis, AnalysisPolicy, Code, Coverage, Diagnostic, Diagnostics, Findings, Format,
    Provenance, Rendered, Status, Structure, StructureView,
};

// The policies and identifiers a call takes.
pub use crate::{
    AltlocPolicy, AmbiguousResidueBoundaryPolicy, AssemblyChoice, AtomRef, AtomSelection, ChainRef,
    ChainSequenceExt, Element, ExecutionContext, Limits, MissingElementPolicy, MissingPolicy,
    ModelChoice, ModelIndex, ModelRef, Namespace, ParseMode, ReadOptions, ResidueRef,
};

// The verbs.
pub use crate::{
    read, read_bytes, read_with_diagnostics, read_with_options, write, write_with_options,
};

#[cfg(feature = "mmcif")]
pub use crate::{WriteOptions, read_document, write_mmcif, write_mmcif_with_options};

#[cfg(feature = "bcif")]
pub use crate::{write_bcif, write_bcif_with_options};

#[cfg(feature = "pdb")]
pub use crate::{PdbHeadersExt, write_pdb};

#[cfg(feature = "chem")]
pub use crate::read_component_dictionary;

#[cfg(feature = "geom")]
pub use crate::{Rigid, transform};

// The bounded batch reader needs a format crate to read with, so it appears
// under exactly the features that give `StructureBatchReader` a variant.
#[cfg(any(
    feature = "mmcif",
    feature = "pdb",
    feature = "bcif",
    feature = "modelcif"
))]
pub use crate::{StructureBatchReader, open_structure_batches};

// Selection.
#[cfg(feature = "query")]
pub use crate::{Groups, Query, QueryStructure};

// The typed operation vocabulary.
#[cfg(all(feature = "analysis", feature = "geom"))]
pub use crate::{
    ContactsRequest, CoordinateInput, CoordinateSlot, ExecutionPlanError, FloatInput, FrameInput,
    GeometryRequest, GeometryValue, IndexInput, PhysicalRequest, PhysicalValue, Plan, PlanInput,
    PlanOperation, PlanResult, PlanResultEntry, PlanValue, RmsdRequest, ScalarInput,
    SelectionRequest, SpatialRequest, SpatialValue, StructureRequest, StructureValue,
};

#[cfg(all(feature = "analysis", feature = "geom", feature = "surface"))]
pub use crate::{MaskInput, SurfaceRequest, SurfaceValue};

#[cfg(all(feature = "analysis", feature = "geom", feature = "traj"))]
pub use crate::{TrajectoryRequest, TrajectoryValue};

// The comparison requests live in `operations`, which needs the two kernels
// its executor is built from.
#[cfg(all(feature = "compare", feature = "analysis", feature = "geom"))]
pub use crate::{ComparisonMetric, ComparisonRequest, ComparisonResult};

#[cfg(feature = "spatial")]
pub use crate::SpatialBackend;

// Structure extensions, so a handle's own methods are callable.
#[cfg(feature = "xtal")]
pub use crate::{AssemblyExt, NcsExt, SymmetryExt};

#[cfg(feature = "modelcif")]
pub use crate::ModelCifExt;
