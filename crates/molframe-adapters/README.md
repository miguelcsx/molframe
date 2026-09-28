# MolFrame adapters

This crate handles integration problems that sit outside the scientific kernels: transferring topology across representation boundaries and acquiring external resources under explicit integrity constraints.

```mermaid
flowchart LR
    Remote["Remote resource"] --> Fetch["bounded fetch"]
    Fetch --> Hash["SHA-256 verification"]
    Hash --> Bytes["VerifiedDownload"]
```

## Topology transfer

`TopologyBatch` is a compact, representation-neutral snapshot built from contiguous columns and batch-local string identifiers. It can be produced or consumed without requiring another library to adopt MolFrame's internal storage layout.

Import validates column lengths, hierarchy offsets, indices, elements, and connectivity before constructing an immutable native structure. Invalid external state is rejected rather than repaired implicitly.

## Resource acquisition

`fetch_verified` retrieves a caller-supplied URL only when the caller also supplies an expected SHA-256 digest, byte ceiling, timeout, and redirect policy. Bytes are published only after integrity verification succeeds. A digest calculated from received bytes establishes content identity, not publisher authenticity.

The crate does not currently ship an authoritative identifier-to-resource manifest. `fetch_identifier` therefore returns the typed `IdentifierFetchError::ProviderUnavailable` diagnostic for non-empty identifiers rather than guessing an archive URL or treating a local digest as an authoritative checksum. Local path and caller-supplied byte loading remain available through the facade's normal readers.
