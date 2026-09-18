//! Allocation-bounded `MessagePack` serialization helpers.

use crate::codec::EncodedData;
use crate::container::{EncodedBlock, EncodedCategory, EncodedColumn, EncodedFile};
use serde::Serialize;
use serde::ser::{Error as _, SerializeSeq, SerializeStruct, Serializer};
use std::cell::RefCell;
use std::io::{self, Write};

pub(crate) fn compact_size<T>(value: &T) -> Result<usize, rmp_serde::encode::Error>
where
    T: Serialize + ?Sized,
{
    let mut counter = ByteCounter::default();
    rmp_serde::encode::write(&mut counter, value)?;
    Ok(counter.bytes)
}

pub(crate) fn consume_to_vec_named(
    value: EncodedFile,
) -> Result<Vec<u8>, rmp_serde::encode::Error> {
    let mut bytes = Vec::new();
    consume_to_writer_named(&mut bytes, value)?;
    Ok(bytes)
}

pub(crate) fn consume_to_writer_named<W: Write>(
    writer: &mut W,
    value: EncodedFile,
) -> Result<(), rmp_serde::encode::Error> {
    rmp_serde::encode::write_named(writer, &Consuming::new(value))
}

pub(crate) struct Consuming<T>(RefCell<Option<T>>);

impl<T> Consuming<T> {
    pub(crate) const fn new(value: T) -> Self {
        Self(RefCell::new(Some(value)))
    }
}

impl<T: SerializeOwned> Serialize for Consuming<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let value = self
            .0
            .borrow_mut()
            .take()
            .ok_or_else(|| S::Error::custom("encoded value was already serialized"))?;
        value.serialize_owned(serializer)
    }
}

trait SerializeOwned {
    fn serialize_owned<S>(self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer;
}

impl SerializeOwned for EncodedFile {
    fn serialize_owned<S>(self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("EncodedFile", 3)?;
        state.serialize_field("version", &self.version)?;
        state.serialize_field("encoder", &self.encoder)?;
        state.serialize_field("dataBlocks", &ConsumingSequence::new(self.data_blocks))?;
        state.end()
    }
}

impl SerializeOwned for EncodedBlock {
    fn serialize_owned<S>(self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("EncodedBlock", 2)?;
        state.serialize_field("header", &self.header)?;
        state.serialize_field("categories", &ConsumingSequence::new(self.categories))?;
        state.end()
    }
}

impl SerializeOwned for EncodedCategory {
    fn serialize_owned<S>(self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("EncodedCategory", 3)?;
        state.serialize_field("name", &self.name)?;
        state.serialize_field("rowCount", &self.row_count)?;
        state.serialize_field("columns", &ConsumingSequence::new(self.columns))?;
        state.end()
    }
}

impl SerializeOwned for EncodedColumn {
    fn serialize_owned<S>(self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let field_count = if self.mask.is_some() { 3 } else { 2 };
        let mut state = serializer.serialize_struct("EncodedColumn", field_count)?;
        state.serialize_field("name", &self.name)?;
        state.serialize_field("data", &Consuming::new(self.data))?;
        if let Some(mask) = self.mask {
            state.serialize_field("mask", &Consuming::new(mask))?;
        }
        state.end()
    }
}

impl SerializeOwned for EncodedData {
    fn serialize_owned<S>(self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("EncodedData", 2)?;
        state.serialize_field("encoding", &self.encoding)?;
        state.serialize_field("data", serde_bytes::Bytes::new(&self.data))?;
        state.end()
    }
}

struct ConsumingSequence<T>(RefCell<Option<Vec<T>>>);

impl<T> ConsumingSequence<T> {
    const fn new(values: Vec<T>) -> Self {
        Self(RefCell::new(Some(values)))
    }
}

impl<T: SerializeOwned> Serialize for ConsumingSequence<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let values = self
            .0
            .borrow_mut()
            .take()
            .ok_or_else(|| S::Error::custom("encoded sequence was already serialized"))?;
        let mut sequence = serializer.serialize_seq(Some(values.len()))?;
        for value in values {
            sequence.serialize_element(&Consuming::new(value))?;
        }
        sequence.end()
    }
}

#[derive(Default)]
struct ByteCounter {
    bytes: usize,
}

impl Write for ByteCounter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        let Some(bytes) = self.bytes.checked_add(buffer.len()) else {
            return Err(io::Error::other("serialized size exceeds usize"));
        };
        self.bytes = bytes;
        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
#[path = "messagepack_tests.rs"]
mod tests;
