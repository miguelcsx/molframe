//! Declarative command-line argument groups.

mod advanced;
mod interactions;
mod mapping;
mod policy;
mod sequence;
mod validate;
mod workflows;

pub(crate) use advanced::{AuditArguments, BatchCommand, FxCommand};
pub(crate) use interactions::{
    HydrogenBondArguments, InteractionCommand, InteractionInput, PlaneFitArguments,
};
pub(crate) use mapping::MappingArguments;
pub(crate) use policy::{AltlocArgument, AssemblyArgument, ModelArgument, NamespaceArgument};
pub(crate) use sequence::{
    KmerOperation, MatrixChoice, PairwiseMode, SequenceCommand, SequenceFormat, TreeMethod,
};
pub(crate) use validate::ValidateArguments;

pub(crate) use workflows::{
    CcdArguments, EnsembleCommand, GeometryCommand, SurfaceArguments, SystemCommand,
    TrajectoryCommand,
};
