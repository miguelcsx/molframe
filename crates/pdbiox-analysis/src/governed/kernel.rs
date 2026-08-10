//! Typed kernel outputs and a closure adapter.

use super::AnalysisDescriptor;
use pdbiox_core::contract::{AnalysisPolicy, Assumption, Coverage, Status};
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::structure::Structure;
use std::marker::PhantomData;

/// One frame's value plus its honest scientific completeness metadata.
#[derive(Clone, Debug)]
pub struct FrameKernelResult<T> {
    /// Kernel value for this frame.
    pub value: T,
    /// Whether the value is defensible under the supplied inputs.
    pub status: Status,
    /// Intended, used, missing and ambiguous inputs.
    pub coverage: Coverage,
    /// Deterministically ordered findings from the kernel.
    pub warnings: Vec<Diagnostic>,
    /// Decisions made while executing this frame.
    pub assumptions: Vec<Assumption>,
}

impl<T> FrameKernelResult<T> {
    /// A frame for which every intended input was used unambiguously.
    #[must_use]
    pub fn complete(value: T, intended: u32) -> Self {
        Self {
            value,
            status: Status::Complete,
            coverage: Coverage::complete(intended),
            warnings: Vec::new(),
            assumptions: Vec::new(),
        }
    }

    /// Builds a result with explicit completeness metadata.
    #[must_use]
    pub fn governed(value: T, status: Status, coverage: Coverage) -> Self {
        Self {
            value,
            status,
            coverage,
            warnings: Vec::new(),
            assumptions: Vec::new(),
        }
    }

    /// Attaches one deterministically ordered finding.
    #[must_use]
    pub fn with_warning(mut self, warning: Diagnostic) -> Self {
        self.warnings.push(warning);
        self
    }

    /// Records one decision made while evaluating this frame.
    #[must_use]
    pub fn with_assumption(mut self, assumption: Assumption) -> Self {
        self.assumptions.push(assumption);
        self
    }
}

/// A Rust kernel that can run on one immutable structure frame.
pub trait StructureKernel: Sync {
    /// Value produced for one frame.
    type Output: Send;
    /// Typed kernel failure.
    type Error: Send;

    /// Stable identity and complete parameter set.
    fn descriptor(&self) -> &AnalysisDescriptor;
    /// Executes exactly one frame under the recorded policy without
    /// reimplementing its scientific logic.
    ///
    /// # Errors
    ///
    /// Returns the original typed kernel error without converting it to a
    /// trajectory or string error.
    fn analyse(
        &self,
        structure: &Structure,
        policy: &AnalysisPolicy,
    ) -> Result<FrameKernelResult<Self::Output>, Self::Error>;

    /// Executes on a policy-materialized frame with its target-to-source atom map.
    ///
    /// Kernels with atom-aligned side inputs override this to project those
    /// inputs without changing their scientific implementation.
    ///
    /// # Errors
    ///
    /// Returns the kernel's typed scientific error.
    fn analyse_mapped(
        &self,
        structure: &Structure,
        policy: &AnalysisPolicy,
        _source_atoms: &[usize],
    ) -> Result<FrameKernelResult<Self::Output>, Self::Error> {
        self.analyse(structure, policy)
    }
}

/// Source-map-aware closure kernel created by [`mapped_structure_kernel`].
#[derive(Debug)]
pub struct MappedClosureStructureKernel<F, T, E> {
    descriptor: AnalysisDescriptor,
    run: F,
    marker: PhantomData<fn() -> (T, E)>,
}

/// Adapts a kernel with atom-aligned side inputs to policy materialization.
#[must_use]
pub fn mapped_structure_kernel<F, T, E>(
    descriptor: AnalysisDescriptor,
    run: F,
) -> MappedClosureStructureKernel<F, T, E>
where
    F: Fn(&Structure, &AnalysisPolicy, &[usize]) -> Result<FrameKernelResult<T>, E> + Sync,
    T: Send,
    E: Send,
{
    MappedClosureStructureKernel {
        descriptor,
        run,
        marker: PhantomData,
    }
}

/// Closure-backed structure kernel created by [`structure_kernel`].
#[derive(Debug)]
pub struct ClosureStructureKernel<F, T, E> {
    descriptor: AnalysisDescriptor,
    run: F,
    marker: PhantomData<fn() -> (T, E)>,
}

/// Adapts an existing Rust kernel without duplicating its implementation.
#[must_use]
pub fn structure_kernel<F, T, E>(
    descriptor: AnalysisDescriptor,
    run: F,
) -> ClosureStructureKernel<F, T, E>
where
    F: Fn(&Structure, &AnalysisPolicy) -> Result<FrameKernelResult<T>, E> + Sync,
    T: Send,
    E: Send,
{
    ClosureStructureKernel {
        descriptor,
        run,
        marker: PhantomData,
    }
}

impl<F, T, E> StructureKernel for ClosureStructureKernel<F, T, E>
where
    F: Fn(&Structure, &AnalysisPolicy) -> Result<FrameKernelResult<T>, E> + Sync,
    T: Send,
    E: Send,
{
    type Output = T;
    type Error = E;

    fn descriptor(&self) -> &AnalysisDescriptor {
        &self.descriptor
    }

    fn analyse(
        &self,
        structure: &Structure,
        policy: &AnalysisPolicy,
    ) -> Result<FrameKernelResult<T>, E> {
        (self.run)(structure, policy)
    }
}

impl<F, T, E> StructureKernel for MappedClosureStructureKernel<F, T, E>
where
    F: Fn(&Structure, &AnalysisPolicy, &[usize]) -> Result<FrameKernelResult<T>, E> + Sync,
    T: Send,
    E: Send,
{
    type Output = T;
    type Error = E;

    fn descriptor(&self) -> &AnalysisDescriptor {
        &self.descriptor
    }

    fn analyse(
        &self,
        structure: &Structure,
        policy: &AnalysisPolicy,
    ) -> Result<FrameKernelResult<T>, E> {
        let identity: Vec<usize> = (0..structure.atom_count() as usize).collect();
        (self.run)(structure, policy, &identity)
    }

    fn analyse_mapped(
        &self,
        structure: &Structure,
        policy: &AnalysisPolicy,
        source_atoms: &[usize],
    ) -> Result<FrameKernelResult<T>, E> {
        (self.run)(structure, policy, source_atoms)
    }
}
