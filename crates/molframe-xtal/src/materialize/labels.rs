// Chain-label ledger that keeps generated labels unique across whole and chunked expansion.
//
// Complexity: amortised O(1) time per label assigned; O(labels) space.

/// Chain-label bookkeeping shared by whole and chunked expansion.
///
/// Labels are tracked as text so a ledger carried from one chunk to the next
/// stays valid even though every chunk owns its own dictionary.
#[derive(Default)]
#[derive(Debug)]
struct LabelLedger {
    used: std::collections::HashSet<Box<str>>,
    counts: IdentityHashMap<SymbolId, u32>,
}

impl LabelLedger {
    /// Generates a stable collision-free label and returns its dictionary symbol.
    ///
    /// The first copy of a source chain keeps the source symbol; later copies
    /// are suffixed `-2`, `-3`, ... and interned only when a collision requires.
    fn assign(
        &mut self,
        dictionary: &mut molframe_core::symbol::Interner,
        source_symbol: SymbolId,
        source: &str,
    ) -> Result<SymbolId, Diagnostic> {
        let mut count = self.counts.get(&source_symbol).copied().unwrap_or(0);
        loop {
            count = count
                .checked_add(1)
                .ok_or_else(|| capacity("chain identifiers"))?;
            let (symbol, text) = if count == 1 {
                (Some(source_symbol), source.to_owned())
            } else {
                (None, format!("{source}-{count}"))
            };
            if self.used.contains(text.as_str()) {
                continue;
            }
            let symbol = match symbol {
                Some(symbol) => symbol,
                None => dictionary
                    .intern(&text)
                    .map_err(|_| capacity("dictionary entries"))?,
            };
            self.used.insert(text.into_boxed_str());
            self.counts.insert(source_symbol, count);
            return Ok(symbol);
        }
    }
}
