//! Stable command-line vocabularies mapped onto public Rust types.

use clap::ValueEnum;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum MetricChoice {
    Lddt,
    TmScore,
    GdtTs,
    GdtHa,
    #[value(name = "dockq")]
    DockQ,
    Qs,
    Cad,
    Ce,
    ContactSimilarity,
    LigandRmsd,
    InterfaceRmsd,
    PocketRmsd,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum EmptyLddtChoice {
    Perfect,
    Error,
}

impl From<EmptyLddtChoice> for molframe::compare::EmptyLddtPolicy {
    fn from(value: EmptyLddtChoice) -> Self {
        match value {
            EmptyLddtChoice::Perfect => Self::Perfect,
            EmptyLddtChoice::Error => Self::Error,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum CompletionShell {
    Bash,
    Elvish,
    Fish,
    PowerShell,
    Zsh,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum ValidationChoice {
    Core,
    Quality,
    Geometry,
    Clashes,
    Completeness,
    Altloc,
    CcdCompleteness,
    Bfactor,
    Valence,
    Ligand,
    Chirality,
    CisPeptide,
    Nucleic,
    Ramachandran,
    Rotamer,
    PlaneRestraints,
    Tls,
    Protein,
    GeometryZ,
    Clashscore,
    NucleicConformation,
}
