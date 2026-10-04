//! Catalogued operations of the `sequence` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[
    Capability {
        name: "align",
        domain: "sequence",
        feature: "sequence",
        inputs: "sequence, sequence",
        result: "Alignment",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "parse_fasta",
        domain: "sequence",
        feature: "sequence",
        inputs: "text",
        result: "list",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "write_fasta",
        domain: "sequence",
        feature: "sequence",
        inputs: "records",
        result: "text",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "kmer_counts",
        domain: "sequence",
        feature: "sequence",
        inputs: "sequence",
        result: "list",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
];
