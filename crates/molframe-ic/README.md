# molframe-ic

Internal-coordinate representation and deterministic Cartesian rebuilding.

This crate converts molecular geometry between Cartesian coordinates and a bond-angle-torsion representation while preserving enough Cartesian seeds to orient each disconnected component.

```mermaid
flowchart LR
    Structure["Structure + bond graph"] --> Forest["Deterministic traversal forest"]
    Forest --> Seeds["Cartesian seeds"]
    Forest --> BAT["Bond · angle · torsion rows"]

    Seeds --> IC["InternalCoordinates"]
    BAT --> IC

    IC --> Measure["measure_bat(frame)"]
    IC --> Rebuild["rebuild_bat(frame)"]
    Rebuild --> XYZ["Cartesian coordinates"]
```

Shallow atoms and disconnected components remain Cartesian seeds. Deeper atoms are represented by references to three previously placed atoms plus bond length, bond angle, and torsion.

The topology of the internal-coordinate forest can be reused across trajectory frames, allowing new BAT values to be measured without repeating graph traversal.

Degenerate or missing geometry is reported rather than reconstructed from invented coordinates.
