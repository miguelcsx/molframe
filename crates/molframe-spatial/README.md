# molframe-spatial

Spatial planning, indexing, and fixed-radius queries.

Rather than exposing one mandatory spatial data structure, `molframe-spatial` separates **what a query needs** from **which backend should answer it**.

```mermaid
flowchart LR
    Request["Selection + cutoff + periodicity"] --> Planner["Spatial planner"]
    Planner --> Backend{"Backend"}
    Backend --> Brute["Brute / SIMD"]
    Backend --> Cell["Cell list"]
    Backend --> KD["k-d tree"]
    Backend --> NL["Neighbor list"]

    Brute --> Query["Fixed-radius operations"]
    Cell --> Query
    KD --> Query
    NL --> Query

    Query --> Pairs["pairs · counts · indices · reductions"]
```

`StructureSpatial` binds spatial execution to one immutable structure snapshot. Compatible indexes can be retained and reused, while cache state remains bounded so query-controlled workloads cannot grow memory indefinitely.

Periodic geometry is explicit through `PeriodicBox` rather than being inferred by individual analyses.

Higher-level domains such as interactions, surfaces, validation, trajectories, and structural comparison reuse this layer instead of implementing independent all-pairs searches.
