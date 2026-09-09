use pdbiox_core::Diagnostic;
use std::fmt;

/// Configuration or storage failure while projecting `ModelCIF` metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ModelCifError {
    /// The configured projection-memory ceiling is outside the supported range.
    InvalidMemoryLimit {
        /// Requested ceiling in bytes.
        limit: usize,
    },
    /// The compact representation cannot fit inside the configured ceiling.
    MemoryLimit {
        /// Peak projection-owned bytes required, including widening buffers.
        required: usize,
        /// Configured ceiling in bytes.
        limit: usize,
    },
    /// A row or dictionary index exceeds the compact public index domain.
    Capacity,
    /// The allocator refused a bounded compact-column allocation.
    Allocation,
    /// Atom-aligned annotations cannot represent different values per dense model.
    MultipleStructureModels {
        /// Number of coordinate models sharing the topology.
        count: usize,
    },
    /// Local pLDDT rows refer to more than one `ModelCIF` model.
    MultiplePlddtModels {
        /// First model encountered in source order.
        first: Box<str>,
        /// Conflicting model encountered later.
        second: Box<str>,
    },
    /// More than one pLDDT row addresses the same residue.
    DuplicatePlddtResidue {
        /// Label chain identifier.
        chain: Box<str>,
        /// Label sequence identifier.
        sequence: i64,
    },
}

impl fmt::Display for ModelCifError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMemoryLimit { limit } => write!(
                formatter,
                "ModelCIF memory limit {limit} must be at least one byte"
            ),
            Self::MemoryLimit { required, limit } => write!(
                formatter,
                "compact ModelCIF metadata needs {required} bytes, exceeding limit {limit}"
            ),
            Self::Capacity => formatter.write_str("ModelCIF row or dictionary capacity exceeded"),
            Self::Allocation => {
                formatter.write_str("compact ModelCIF column allocation was refused")
            }
            Self::MultipleStructureModels { count } => write!(
                formatter,
                "cannot attach one pLDDT column to {count} coordinate models"
            ),
            Self::MultiplePlddtModels { first, second } => write!(
                formatter,
                "pLDDT rows address multiple ModelCIF models: {first} and {second}"
            ),
            Self::DuplicatePlddtResidue { chain, sequence } => {
                write!(formatter, "duplicate pLDDT for residue {chain}:{sequence}")
            }
        }
    }
}

impl std::error::Error for ModelCifError {}

/// Syntax or resource failure while reading compact `ModelCIF` directly.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ModelCifReadError {
    /// CIF syntax findings prevented a complete projection.
    Syntax(Vec<Diagnostic>),
    /// The compact projection violated its capacity or memory policy.
    Projection(ModelCifError),
}

impl fmt::Display for ModelCifReadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Syntax(findings) => {
                write!(formatter, "CIF syntax prevented ModelCIF projection")?;
                if let Some(finding) = findings.first() {
                    write!(formatter, ": {finding}")?;
                }
                Ok(())
            }
            Self::Projection(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for ModelCifReadError {}

/// Resource policy for `ModelCIF` projection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModelCifOptions {
    /// Maximum projection-owned bytes for compact `ma_*` storage.
    pub memory_limit_bytes: usize,
}

impl ModelCifOptions {
    /// Default peak projection-memory ceiling.
    pub const DEFAULT_MEMORY_LIMIT_BYTES: usize = 100_000_000;

    /// Creates the default bounded projection policy.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            memory_limit_bytes: Self::DEFAULT_MEMORY_LIMIT_BYTES,
        }
    }

    /// Sets the peak projection-memory ceiling.
    #[must_use]
    pub const fn with_memory_limit(mut self, bytes: usize) -> Self {
        self.memory_limit_bytes = bytes;
        self
    }

    pub(crate) fn validate(self) -> Result<Self, ModelCifError> {
        if self.memory_limit_bytes == 0 {
            return Err(ModelCifError::InvalidMemoryLimit {
                limit: self.memory_limit_bytes,
            });
        }
        Ok(self)
    }
}

impl Default for ModelCifOptions {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct MemoryBudget {
    used: usize,
    limit: usize,
}

impl MemoryBudget {
    pub(crate) const fn new(limit: usize) -> Self {
        Self { used: 0, limit }
    }

    pub(crate) fn claim(&mut self, bytes: usize) -> Result<(), ModelCifError> {
        let required = self
            .used
            .checked_add(bytes)
            .ok_or(ModelCifError::Capacity)?;
        if required > self.limit {
            return Err(ModelCifError::MemoryLimit {
                required,
                limit: self.limit,
            });
        }
        self.used = required;
        Ok(())
    }
}
