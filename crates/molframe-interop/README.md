# molframe-interop

Zero-copy and structured interoperability for MolFrame data.

`molframe-interop` connects MolFrame's immutable structural model to external data ecosystems without introducing a second molecular representation.

It provides columnar tables, streaming interfaces, tensor exchange, graph projections, and dataset utilities built directly from `molframe-core::Structure`.

```mermaid
flowchart LR
    Structure["MolFrame Structure"]

    Structure --> Arrow["Arrow tables"]
    Structure --> Graph["Graph projection"]
    Structure --> DLPack["DLPack tensor"]
    Structure --> Dataset["Dataset manifests"]

    Arrow --> Stream["Arrow C Stream"]
    Arrow --> IPC["Arrow IPC"]
    Arrow --> Parquet["Parquet"]

    Graph --> Tensor["Graph / tensor frameworks"]
    DLPack --> Tensor

    Dataset --> Filter["Filter"]
    Dataset --> Split["Deterministic split"]
```

## Architecture

`Structure` remains the source of truth.

The crate does not copy structural semantics into framework-specific object models. Instead, it projects the same immutable snapshot into representations designed for analytics, tensor libraries, graph frameworks, or persistent columnar storage.

### Columnar interoperability

`AtomTable`, `ResidueTable`, `ChainTable`, and `BondTable` expose structural data through Apache Arrow.

Atom record batches follow MolFrame's internal chunk boundaries:

```mermaid
flowchart LR
    Chunks["Structure atom chunks"]
        --> Batch1["Arrow RecordBatch"]
        --> Consumer["Arrow consumer"]

    Chunks --> Batch2["Arrow RecordBatch"]
```

This allows consumers to process large structures incrementally rather than requiring one fully materialized table.

`ArrowStream` exposes the same model through the Arrow C Stream Interface and materializes at most one batch as the consumer pulls data.

Arrow-backed tables can also be written as:

- Arrow IPC
- Parquet

Schema metadata can retain MolFrame-specific information alongside the tabular representation.

## Tensor interoperability

Coordinate arrays can be exported through DLPack for tensor consumers.

Unlike Arrow, DLPack does not provide a read-only ownership contract: a consumer may legally mutate the tensor it receives. MolFrame structures are immutable snapshots, so exposing their coordinate storage directly would violate the structural model.

For that reason, DLPack export creates independently owned contiguous storage:

```mermaid
flowchart LR
    Structure["Immutable coordinates"]
        --> Copy["Owned contiguous copy"]
        --> DLPack["DLManagedTensor"]
        --> Consumer["Tensor consumer"]
```

The copy is therefore an intentional ownership boundary rather than an implementation limitation.

## Graph projection

Molecular structures can be projected into deterministic graph representations suitable for graph-oriented consumers.

Graph construction derives nodes, edges, and features from the existing structural model instead of maintaining an independent graph-backed molecule representation.

```mermaid
flowchart LR
    Structure --> Nodes["Nodes"]
    Structure --> Edges["Edges"]
    Structure --> Features["Features"]

    Nodes --> Graph["Graph representation"]
    Edges --> Graph
    Features --> Graph
```

This layer defines **data interchange**, not graph-learning algorithms.

## Dataset layer

Manifest-backed datasets provide lazy organization of collections of structures without loading the complete collection into memory.

Dataset utilities include filtering and deterministic splitting for reproducible downstream workflows.

```mermaid
flowchart LR
    Manifest["Dataset manifest"]
        --> Dataset["Lazy Dataset"]

    Dataset --> Filtered["Filtered view"]
    Dataset --> Train["Train"]
    Dataset --> Validation["Validation"]
    Dataset --> Test["Test"]
```

Dataset management remains independent of any particular machine-learning framework.

## Boundary

This crate owns representations used to move MolFrame data across ecosystem boundaries.

It does **not** implement:

- machine-learning models;
- training or inference;
- tensor operations;
- graph neural networks;
- molecular structure storage;
- scientific analysis kernels.

Those consumers operate on data projected from the shared MolFrame `Structure`.
