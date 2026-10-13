//! Shared accounting for bounded execution memory.

use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Once, OnceLock};

/// Smallest automatic execution budget, and the legacy fixed default: 100 MB.
///
/// The automatic default ([`MemoryBudget::auto`]) never falls below this. A
/// caller who provisioned the machine is the only party who knows what the
/// working set may be, so a budget larger or smaller than this is accepted
/// without argument.
pub const DEFAULT_MEMORY_BUDGET_BYTES: usize = 100_000_000;

/// Environment variable that overrides the automatic budget (`4GiB`, `512MB`, bytes).
pub const MEMORY_BUDGET_ENV: &str = "MOLFRAME_MEMORY_BUDGET";

/// Largest automatic budget: 512 GiB. Explicit budgets are not capped.
pub const MAX_AUTOMATIC_BUDGET_BYTES: usize = 512 << 30;

/// Budget assumed when the platform cannot report its memory (1 GiB).
const UNKNOWN_PLATFORM_BYTES: usize = 1 << 30;

/// The text every memory refusal appends so the caller can act on it.
pub const MEMORY_BUDGET_REMEDY: &str = "raise the memory budget: Python \
ExecutionContext(memory_budget=...), Rust ExecutionContext::builder().memory_budget(MemoryBudget::new(..)?), \
CLI --memory-budget, or environment variable MOLFRAME_MEMORY_BUDGET";

/// A validated upper bound for memory retained by one execution context.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MemoryBudget {
    bytes: usize,
}

impl MemoryBudget {
    /// Creates a non-zero budget.
    ///
    /// There is no compile-time maximum. A ceiling written into a library is a
    /// ceiling nobody can raise, and a caller reading a hundred-gigabyte entry
    /// needs a working set no constant here could have anticipated.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Zero`] when `bytes` is zero, because a
    /// zero-byte budget cannot make forward progress.
    pub const fn new(bytes: usize) -> Result<Self, MemoryBudgetError> {
        if bytes == 0 {
            return Err(MemoryBudgetError::Zero);
        }
        Ok(Self { bytes })
    }

    /// Returns the maximum number of retained bytes.
    #[must_use]
    pub const fn bytes(self) -> usize {
        self.bytes
    }

    /// The budget this process uses when the caller sets none.
    ///
    /// `MOLFRAME_MEMORY_BUDGET` wins when it parses (see [`MemoryBudget::parse`]);
    /// an unparsable value is reported once on stderr and ignored. Otherwise the
    /// budget is [`MemoryBudget::for_physical_memory`] of the physical memory this
    /// process may use: the container limit when one is set, available memory on
    /// Linux, total memory on macOS. A platform that cannot report memory gets 1 GiB.
    #[must_use]
    pub fn auto() -> Self {
        static DETECTED: OnceLock<MemoryBudget> = OnceLock::new();
        if let Ok(text) = std::env::var(MEMORY_BUDGET_ENV) {
            match Self::parse(&text) {
                Ok(budget) => return budget,
                Err(_) => warn_invalid_env(&text),
            }
        }
        *DETECTED.get_or_init(|| match detect_physical_bytes() {
            Some(bytes) => Self::for_physical_memory(bytes),
            None => Self {
                bytes: UNKNOWN_PLATFORM_BYTES,
            },
        })
    }

    /// Forty percent of `physical_bytes`, clamped to 100 MB through 512 GiB.
    #[must_use]
    pub fn for_physical_memory(physical_bytes: usize) -> Self {
        let share = physical_bytes / 5 * 2;
        Self {
            bytes: share.clamp(DEFAULT_MEMORY_BUDGET_BYTES, MAX_AUTOMATIC_BUDGET_BYTES),
        }
    }

    /// Parses a byte count such as `123456`, `512MB`, `4GiB` or `1.5G`.
    ///
    /// `KB`, `MB`, `GB`, `TB` are powers of 1000; `K`, `KiB` and so on, and a bare
    /// `K`/`M`/`G`/`T`, are powers of 1024.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Zero`] for zero and
    /// [`MemoryBudgetError::Unparsable`] for anything else that is not a size.
    pub fn parse(text: &str) -> Result<Self, MemoryBudgetError> {
        let trimmed = text.trim();
        let split = match trimmed.find(|c: char| !(c.is_ascii_digit() || c == '.')) {
            Some(split) => split,
            None => trimmed.len(),
        };
        let (number, unit) = trimmed.split_at(split);
        let multiplier: f64 = match unit.trim().to_ascii_lowercase().as_str() {
            "" | "b" => 1.0,
            "k" | "kib" => 1024.0,
            "m" | "mib" => 1024.0 * 1024.0,
            "g" | "gib" => 1024.0 * 1024.0 * 1024.0,
            "t" | "tib" => 1024.0_f64.powi(4),
            "kb" => 1e3,
            "mb" => 1e6,
            "gb" => 1e9,
            "tb" => 1e12,
            _ => return Err(MemoryBudgetError::Unparsable),
        };
        let value: f64 = number.parse().map_err(|_| MemoryBudgetError::Unparsable)?;
        let bytes = value * multiplier;
        #[allow(clippy::cast_precision_loss)]
        let ceiling = usize::MAX as f64;
        if !bytes.is_finite() || bytes >= ceiling {
            return Err(MemoryBudgetError::Unparsable);
        }
        // Truncation is the intent: a fractional byte is not a byte.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        Self::new(bytes as usize)
    }
}

impl Default for MemoryBudget {
    /// The machine-scaled budget of [`MemoryBudget::auto`].
    fn default() -> Self {
        Self::auto()
    }
}

fn warn_invalid_env(text: &str) {
    static WARNED: Once = Once::new();
    WARNED.call_once(|| {
        eprintln!(
            "molframe: ignoring {MEMORY_BUDGET_ENV}={text:?}: expected a positive size such as 4GiB or 512MB"
        );
    });
}

/// Physical bytes this process may use, when the platform says so cheaply.
fn detect_physical_bytes() -> Option<usize> {
    let limit = cgroup_limit();
    let physical = if cfg!(target_os = "linux") {
        std::fs::read_to_string("/proc/meminfo")
            .ok()
            .and_then(|text| {
                meminfo_kib(&text, "MemAvailable:").or_else(|| meminfo_kib(&text, "MemTotal:"))
            })
            .map(|kib| kib.saturating_mul(1024))
    } else if cfg!(target_os = "macos") {
        std::process::Command::new("sysctl")
            .args(["-n", "hw.memsize"])
            .output()
            .ok()
            .and_then(|out| String::from_utf8(out.stdout).ok())
            .and_then(|text| text.trim().parse::<usize>().ok())
    } else {
        None
    };
    match (physical, limit) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (a, b) => a.or(b),
    }
}

fn meminfo_kib(text: &str, key: &str) -> Option<usize> {
    text.lines()
        .find_map(|line| line.strip_prefix(key))
        .and_then(|rest| rest.split_whitespace().next())
        .and_then(|value| value.parse().ok())
}

/// The container memory limit (cgroup v2, then v1), when one is set.
fn cgroup_limit() -> Option<usize> {
    [
        "/sys/fs/cgroup/memory.max",
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
    ]
    .iter()
    .find_map(|path| std::fs::read_to_string(path).ok())
    .and_then(|text| text.trim().parse::<u64>().ok())
    // cgroup v1 reports "no limit" as a number near 2^63. The comparison is
    // done in u64 so a 32-bit target does not overflow the shift.
    .filter(|limit| *limit < (1_u64 << 60))
    .and_then(|limit| usize::try_from(limit).ok())
}

/// Failure to construct or reserve an execution-memory budget.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MemoryBudgetError {
    /// A zero-byte budget cannot make forward progress.
    Zero,
    /// The text is not a positive size such as `4GiB`.
    Unparsable,
    /// The context has insufficient unreserved capacity.
    Exhausted {
        /// Bytes requested by the operation.
        requested: usize,
        /// Bytes still available at the time of the attempt.
        available: usize,
    },
}

impl fmt::Display for MemoryBudgetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zero => formatter.write_str("execution memory budget must be non-zero"),
            Self::Unparsable => {
                formatter.write_str("memory budget must be a positive size such as 4GiB or 512MB")
            }
            Self::Exhausted {
                requested,
                available,
            } => write!(
                formatter,
                "execution requested {requested} bytes with {available} bytes available; \
                 {MEMORY_BUDGET_REMEDY}"
            ),
        }
    }
}

impl std::error::Error for MemoryBudgetError {}

#[derive(Debug)]
pub(crate) struct MemoryAccount {
    budget: MemoryBudget,
    reserved: AtomicUsize,
    peak_reserved: AtomicUsize,
    live_batches: AtomicUsize,
    peak_live_batches: AtomicUsize,
}

impl MemoryAccount {
    pub(crate) fn new(budget: MemoryBudget) -> Self {
        Self {
            budget,
            reserved: AtomicUsize::new(0),
            peak_reserved: AtomicUsize::new(0),
            live_batches: AtomicUsize::new(0),
            peak_live_batches: AtomicUsize::new(0),
        }
    }

    pub(crate) const fn budget(&self) -> MemoryBudget {
        self.budget
    }

    pub(crate) fn reserve(
        account: &Arc<Self>,
        bytes: usize,
    ) -> Result<MemoryReservation, MemoryBudgetError> {
        let mut current = account.reserved.load(Ordering::Acquire);
        loop {
            let Some(next) = current.checked_add(bytes) else {
                return Err(MemoryBudgetError::Exhausted {
                    requested: bytes,
                    available: account.budget.bytes().saturating_sub(current),
                });
            };
            if next > account.budget.bytes() {
                return Err(MemoryBudgetError::Exhausted {
                    requested: bytes,
                    available: account.budget.bytes().saturating_sub(current),
                });
            }
            match account.reserved.compare_exchange_weak(
                current,
                next,
                Ordering::AcqRel,
                Ordering::Acquire,
            ) {
                Ok(_) => {
                    account.peak_reserved.fetch_max(next, Ordering::AcqRel);
                    return Ok(MemoryReservation {
                        account: Arc::clone(account),
                        bytes,
                    });
                }
                Err(observed) => current = observed,
            }
        }
    }

    pub(crate) fn reserved_bytes(&self) -> usize {
        self.reserved.load(Ordering::Acquire)
    }

    pub(crate) fn peak_reserved_bytes(&self) -> usize {
        self.peak_reserved.load(Ordering::Acquire)
    }

    pub(crate) fn live_batches(&self) -> usize {
        self.live_batches.load(Ordering::Acquire)
    }

    pub(crate) fn peak_live_batches(&self) -> usize {
        self.peak_live_batches.load(Ordering::Acquire)
    }

    fn begin_batch(account: &Arc<Self>) -> BatchGuard {
        let live = account.live_batches.fetch_add(1, Ordering::AcqRel) + 1;
        account.peak_live_batches.fetch_max(live, Ordering::AcqRel);
        BatchGuard {
            account: Arc::clone(account),
        }
    }
}

/// Exclusive accounting for bytes retained by an operation.
///
/// Dropping the value returns its capacity to the context. Reservations cannot
/// be cloned, preventing one charge from representing multiple owners.
#[derive(Debug)]
pub struct MemoryReservation {
    account: Arc<MemoryAccount>,
    bytes: usize,
}

impl MemoryReservation {
    /// Reports whether this reservation belongs to the supplied execution account.
    #[must_use]
    pub fn belongs_to(&self, context: &super::ExecutionContext) -> bool {
        Arc::ptr_eq(&self.account, &context.memory)
    }

    /// Returns the charged byte count.
    #[must_use]
    pub const fn bytes(&self) -> usize {
        self.bytes
    }

    /// Reduces this reservation to the exact retained byte count.
    pub fn shrink_to(&mut self, bytes: usize) {
        if bytes >= self.bytes {
            return;
        }
        let released = self.bytes - bytes;
        self.bytes = bytes;
        self.account.reserved.fetch_sub(released, Ordering::AcqRel);
    }

    /// Extends this reservation before its owner grows an allocation.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Exhausted`] without changing this
    /// reservation when the additional bytes do not fit.
    pub fn try_grow(&mut self, additional: usize) -> Result<(), MemoryBudgetError> {
        if additional == 0 {
            return Ok(());
        }
        let mut extra = MemoryAccount::reserve(&self.account, additional)?;
        let Some(bytes) = self.bytes.checked_add(additional) else {
            return Err(MemoryBudgetError::Exhausted {
                requested: additional,
                available: 0,
            });
        };
        self.bytes = bytes;
        extra.bytes = 0;
        Ok(())
    }

    pub(crate) fn begin_batch(&self) -> BatchGuard {
        MemoryAccount::begin_batch(&self.account)
    }
}

impl Drop for MemoryReservation {
    fn drop(&mut self) {
        self.account
            .reserved
            .fetch_sub(self.bytes, Ordering::AcqRel);
    }
}

/// Exactly-once accounting for one live batch lease.
#[derive(Debug)]
pub(crate) struct BatchGuard {
    account: Arc<MemoryAccount>,
}

impl Drop for BatchGuard {
    fn drop(&mut self) {
        self.account.live_batches.fetch_sub(1, Ordering::AcqRel);
    }
}

#[cfg(test)]
#[path = "memory_tests.rs"]
mod tests;

crate::diagnostic_from!(
    MemoryBudgetError,
    |_error| crate::Code::E7001,
    |diagnostic, error| match error {
        MemoryBudgetError::Zero => diagnostic.with_context("budget", "zero"),
        MemoryBudgetError::Unparsable => diagnostic.with_context("budget", "unparsable"),
        MemoryBudgetError::Exhausted {
            requested,
            available,
        } => diagnostic
            .with_context("requested_bytes", requested.to_string())
            .with_context("available_bytes", available.to_string()),
    }
);
