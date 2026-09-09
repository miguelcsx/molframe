use super::packed::{GrowingIndices, GrowingIntegers, push, reserved};
use crate::{
    CompactColumn, DictionaryColumn, DictionaryEntry, IntegerColumn, ModelCifError, NumberColumn,
};
use pdbiox_cif::{CifScalar, CifValue};
use std::collections::HashMap;
use std::mem::size_of;
use std::sync::Arc;

#[derive(Clone, Copy)]
pub(super) struct PushMemory {
    pub(super) retained: usize,
    pub(super) peak: usize,
}

pub(super) enum ColumnBuilder {
    Dictionary(DictionaryBuilder),
    Integer(IntegerBuilder),
    Number(NumberBuilder),
}

impl ColumnBuilder {
    pub(super) fn new(dictionary: bool) -> Self {
        if dictionary {
            Self::Dictionary(DictionaryBuilder::new())
        } else {
            Self::Integer(IntegerBuilder::new())
        }
    }

    pub(super) fn len(&self) -> usize {
        match self {
            Self::Dictionary(values) => values.len(),
            Self::Integer(values) => values.len(),
            Self::Number(values) => values.len(),
        }
    }

    pub(super) fn live_bytes(&self) -> usize {
        match self {
            Self::Dictionary(values) => values.live_bytes(),
            Self::Integer(values) => values.live_bytes(),
            Self::Number(values) => values.live_bytes(),
        }
    }

    pub(super) fn push(&mut self, value: CifScalar<'_>) -> Result<PushMemory, ModelCifError> {
        let before = self.live_bytes();
        let before_buffers = self.buffer_bytes();
        let mut promoted = false;
        if matches!(value, CifScalar::Text(_)) && !matches!(self, Self::Dictionary(_)) {
            self.promote_dictionary()?;
            promoted = true;
        } else if matches!(value, CifScalar::Float(_)) && matches!(self, Self::Integer(_)) {
            self.promote_number()?;
            promoted = true;
        }
        match self {
            Self::Dictionary(values) => values.push(value)?,
            Self::Integer(values) => values.push(value)?,
            Self::Number(values) => values.push(value)?,
        }
        let retained = self.live_bytes();
        let buffers_changed = before_buffers != self.buffer_bytes();
        let peak = if promoted {
            before.saturating_add(retained.saturating_mul(2))
        } else if buffers_changed {
            before.saturating_add(retained)
        } else {
            retained
        };
        Ok(PushMemory { retained, peak })
    }

    pub(super) fn finish(self) -> CompactColumn {
        match self {
            Self::Dictionary(values) => CompactColumn::Dictionary(values.finish()),
            Self::Integer(values) => CompactColumn::Integer(values.finish()),
            Self::Number(values) => CompactColumn::Number(values.finish()),
        }
    }

    fn promote_dictionary(&mut self) -> Result<(), ModelCifError> {
        let source = std::mem::replace(self, Self::Dictionary(DictionaryBuilder::new()));
        let mut target = DictionaryBuilder::new();
        match source {
            Self::Integer(values) => {
                for row in 0..values.len() {
                    target.push(values.scalar(row))?;
                }
            }
            Self::Number(values) => {
                for row in 0..values.len() {
                    target.push(values.scalar(row))?;
                }
            }
            Self::Dictionary(values) => target = values,
        }
        *self = Self::Dictionary(target);
        Ok(())
    }

    fn promote_number(&mut self) -> Result<(), ModelCifError> {
        let source = std::mem::replace(self, Self::Number(NumberBuilder::new()));
        let Self::Integer(values) = source else {
            *self = source;
            return Ok(());
        };
        let mut target = NumberBuilder::with_capacity(values.len())?;
        for row in 0..values.len() {
            target.push(values.scalar(row))?;
        }
        *self = Self::Number(target);
        Ok(())
    }

    fn buffer_bytes(&self) -> usize {
        match self {
            Self::Dictionary(values) => values.buffer_bytes(),
            Self::Integer(values) => values.live_bytes(),
            Self::Number(values) => values.live_bytes(),
        }
    }
}

pub(super) struct DictionaryBuilder {
    entries: Vec<DictionaryEntry>,
    indices: GrowingIndices,
    text: HashMap<Arc<str>, u32>,
    integers: HashMap<i64, u32>,
    floats: HashMap<u64, u32>,
    unknown: Option<u32>,
    inapplicable: Option<u32>,
    last: Option<u32>,
    entry_payload_bytes: usize,
}

impl DictionaryBuilder {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
            indices: GrowingIndices::new(),
            text: HashMap::new(),
            integers: HashMap::new(),
            floats: HashMap::new(),
            unknown: None,
            inapplicable: None,
            last: None,
            entry_payload_bytes: 0,
        }
    }

    fn len(&self) -> usize {
        self.indices.len()
    }

    fn push(&mut self, value: CifScalar<'_>) -> Result<(), ModelCifError> {
        let index = if let Some(index) = self.repeated_index(value) {
            index
        } else {
            match value {
                CifScalar::Text(value) => self.text_index(value)?,
                CifScalar::Integer(value) => self.integer_index(value)?,
                CifScalar::Float(value) => self.float_index(value)?,
                CifScalar::Unknown => self.sentinel_index(false)?,
                CifScalar::Inapplicable => self.sentinel_index(true)?,
            }
        };
        self.last = Some(index);
        self.indices.push(index)
    }

    fn repeated_index(&self, value: CifScalar<'_>) -> Option<u32> {
        let index = self.last?;
        let entry = self.entries.get(usize::try_from(index).ok()?)?;
        let matches = match (&entry.value, value) {
            (CifValue::Inapplicable, CifScalar::Inapplicable)
            | (CifValue::Unknown, CifScalar::Unknown) => true,
            (CifValue::Text(left), CifScalar::Text(right)) => left.as_ref() == right,
            (CifValue::Integer(left), CifScalar::Integer(right)) => *left == right,
            (CifValue::Float(left), CifScalar::Float(right)) => left.to_bits() == right.to_bits(),
            _ => false,
        };
        matches.then_some(index)
    }

    fn finish(self) -> DictionaryColumn {
        DictionaryColumn {
            entries: self.entries,
            indices: self.indices.finish(),
        }
    }

    fn live_bytes(&self) -> usize {
        self.buffer_bytes().saturating_add(self.entry_payload_bytes)
    }

    fn buffer_bytes(&self) -> usize {
        self.entries
            .capacity()
            .saturating_mul(size_of::<DictionaryEntry>())
            .saturating_add(self.indices.live_bytes())
            .saturating_add(table_bytes(&self.text))
            .saturating_add(table_bytes(&self.integers))
            .saturating_add(table_bytes(&self.floats))
    }

    fn text_index(&mut self, value: &str) -> Result<u32, ModelCifError> {
        if let Some(index) = self.text.get(value).copied() {
            return Ok(index);
        }
        self.text
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        let retained: Arc<str> = value.into();
        let entry = CifValue::Text(Arc::clone(&retained));
        let index = self.insert(&entry)?;
        self.text.insert(retained, index);
        Ok(index)
    }

    fn integer_index(&mut self, value: i64) -> Result<u32, ModelCifError> {
        if let Some(index) = self.integers.get(&value).copied() {
            return Ok(index);
        }
        self.integers
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        let index = self.insert(&CifValue::Integer(value))?;
        self.integers.insert(value, index);
        Ok(index)
    }

    fn float_index(&mut self, value: f64) -> Result<u32, ModelCifError> {
        let bits = value.to_bits();
        if let Some(index) = self.floats.get(&bits).copied() {
            return Ok(index);
        }
        self.floats
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        let index = self.insert(&CifValue::Float(value))?;
        self.floats.insert(bits, index);
        Ok(index)
    }

    fn sentinel_index(&mut self, inapplicable: bool) -> Result<u32, ModelCifError> {
        let existing = if inapplicable {
            self.inapplicable
        } else {
            self.unknown
        };
        if let Some(index) = existing {
            return Ok(index);
        }
        let value = if inapplicable {
            CifValue::Inapplicable
        } else {
            CifValue::Unknown
        };
        let index = self.insert(&value)?;
        if inapplicable {
            self.inapplicable = Some(index);
        } else {
            self.unknown = Some(index);
        }
        Ok(index)
    }

    fn insert(&mut self, value: &CifValue) -> Result<u32, ModelCifError> {
        let index = u32::try_from(self.entries.len()).map_err(|_| ModelCifError::Capacity)?;
        let entry = DictionaryEntry::new(value);
        let payload = entry
            .numeric_text
            .as_ref()
            .map_or(0, |value| value.len())
            .saturating_add(entry.value.as_str().map_or(0, str::len));
        self.entry_payload_bytes = self
            .entry_payload_bytes
            .checked_add(payload)
            .ok_or(ModelCifError::Capacity)?;
        push(&mut self.entries, entry)?;
        Ok(index)
    }
}

pub(super) struct IntegerBuilder {
    values: GrowingIntegers,
    states: Option<Vec<u8>>,
}

impl IntegerBuilder {
    const fn new() -> Self {
        Self {
            values: GrowingIntegers::new(),
            states: None,
        }
    }

    fn len(&self) -> usize {
        self.values.len()
    }

    fn live_bytes(&self) -> usize {
        self.values
            .live_bytes()
            .saturating_add(self.states.as_ref().map_or(0, Vec::capacity))
    }

    fn push(&mut self, value: CifScalar<'_>) -> Result<(), ModelCifError> {
        let (number, state) = match value {
            CifScalar::Integer(value) => (value, 0),
            CifScalar::Unknown => (0, 1),
            CifScalar::Inapplicable => (0, 2),
            CifScalar::Text(_) | CifScalar::Float(_) => return Err(ModelCifError::Capacity),
        };
        self.push_state(state)?;
        self.values.push(number)
    }

    fn scalar(&self, row: usize) -> CifScalar<'static> {
        match self.states.as_ref().and_then(|states| states.get(row)) {
            Some(1) => CifScalar::Unknown,
            Some(2) => CifScalar::Inapplicable,
            _ => CifScalar::Integer(match self.values.get(row) {
                Some(value) => value,
                None => 0,
            }),
        }
    }

    fn push_state(&mut self, state: u8) -> Result<(), ModelCifError> {
        if state != 0 && self.states.is_none() {
            let mut states = reserved(self.len())?;
            states.resize(self.len(), 0);
            self.states = Some(states);
        }
        if let Some(states) = &mut self.states {
            push(states, state)?;
        }
        Ok(())
    }

    fn finish(self) -> IntegerColumn {
        IntegerColumn {
            values: self.values.finish(),
            states: self.states,
        }
    }
}

pub(super) struct NumberBuilder {
    payload: Vec<u64>,
    kinds: Option<Vec<u8>>,
}

impl NumberBuilder {
    const fn new() -> Self {
        Self {
            payload: Vec::new(),
            kinds: None,
        }
    }

    fn with_capacity(capacity: usize) -> Result<Self, ModelCifError> {
        Ok(Self {
            payload: reserved(capacity)?,
            kinds: None,
        })
    }

    fn len(&self) -> usize {
        self.payload.len()
    }

    fn live_bytes(&self) -> usize {
        self.payload
            .capacity()
            .saturating_mul(size_of::<u64>())
            .saturating_add(self.kinds.as_ref().map_or(0, Vec::capacity))
    }

    fn push(&mut self, value: CifScalar<'_>) -> Result<(), ModelCifError> {
        let (payload, kind) = match value {
            CifScalar::Float(value) => (value.to_bits(), 0),
            CifScalar::Integer(value) => (u64::from_ne_bytes(value.to_ne_bytes()), 1),
            CifScalar::Unknown => (0, 2),
            CifScalar::Inapplicable => (0, 3),
            CifScalar::Text(_) => return Err(ModelCifError::Capacity),
        };
        self.push_kind(kind)?;
        push(&mut self.payload, payload)
    }

    fn scalar(&self, row: usize) -> CifScalar<'static> {
        let bits = match self.payload.get(row).copied() {
            Some(bits) => bits,
            None => 0,
        };
        match self.kinds.as_ref().and_then(|kinds| kinds.get(row)) {
            Some(1) => CifScalar::Integer(i64::from_ne_bytes(bits.to_ne_bytes())),
            Some(2) => CifScalar::Unknown,
            Some(3) => CifScalar::Inapplicable,
            _ => CifScalar::Float(f64::from_bits(bits)),
        }
    }

    fn push_kind(&mut self, kind: u8) -> Result<(), ModelCifError> {
        if kind != 0 && self.kinds.is_none() {
            let mut kinds = reserved(self.len())?;
            kinds.resize(self.len(), 0);
            self.kinds = Some(kinds);
        }
        if let Some(kinds) = &mut self.kinds {
            push(kinds, kind)?;
        }
        Ok(())
    }

    fn finish(self) -> NumberColumn {
        NumberColumn {
            payload: self.payload,
            kinds: self.kinds,
        }
    }
}

fn table_bytes<K, V>(values: &HashMap<K, V>) -> usize {
    values.capacity().saturating_mul(
        size_of::<(K, V)>()
            .saturating_add(size_of::<usize>())
            .saturating_add(1),
    )
}
