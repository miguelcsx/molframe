//! Declarative command-line argument groups.

mod advanced;
mod policy;
mod sequence;
mod workflows;

pub(crate) use advanced::{AuditArguments, BatchCommand, FxCommand};
pub(crate) use policy::{AltlocArgument, AssemblyArgument, ModelArgument, NamespaceArgument};
pub(crate) use sequence::{
    KmerOperation, MatrixChoice, PairwiseMode, SequenceCommand, SequenceFormat, TreeMethod,
};

pub(crate) use workflows::{
    CcdArguments, EnsembleCommand, GeometryCommand, SurfaceArguments, SystemCommand,
    TrajectoryCommand,
};
