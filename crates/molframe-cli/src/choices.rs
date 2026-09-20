//! Stable command-line vocabularies mapped onto public Rust types.

use clap::ValueEnum;

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum RadiusChoice {
    Bondi,
    AmberUnited,
    Charmm,
    Alvarez,
}

impl From<RadiusChoice> for molframe::chemistry::RadiusSet {
    fn from(value: RadiusChoice) -> Self {
        match value {
            RadiusChoice::Bondi => Self::Bondi,
            RadiusChoice::AmberUnited => Self::AmberUnited,
            RadiusChoice::Charmm => Self::Charmm,
            RadiusChoice::Alvarez => Self::Alvarez,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum MetricChoice {
    Lddt,
    TmScore,
    GdtTs,
    GdtHa,
    #[value(name = "dockq")]
    DockQ,
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
}
