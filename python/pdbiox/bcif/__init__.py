"""Native BinaryCIF codecs, documents, and structure round trips."""

from .._native import (
    BcifReader,
    BinaryDocument,
    DataType,
    Decoded,
    Document,
    EncodedData,
    Encoding,
    decode,
    encode_floats,
    encode_integers,
    encode_interval,
    encode_strings,
    read_bcif as read,
    read_bcif_document as read_document,
    read_bcif_with_document as read_with_document,
    write_bcif as write_structure,
    write_bcif_document as write_document,
    write_bcif_with_options as write_structure_with_options,
)

__all__ = [
    "BcifReader", "BinaryDocument", "DataType", "Decoded", "Document", "EncodedData", "Encoding",
    "decode", "encode_floats", "encode_integers", "encode_interval", "encode_strings",
    "read", "read_document", "read_with_document", "write_document", "write_structure",
    "write_structure_with_options",
]
