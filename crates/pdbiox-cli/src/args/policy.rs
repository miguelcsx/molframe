//! Typed command-line policy overrides with strict, descriptive parsers.

use clap::ValueEnum;
use std::str::FromStr;

#[derive(Clone, Debug)]
pub(crate) enum AssemblyArgument {
    AsymmetricUnit,
    Biological(Box<str>),
    Crystal { radius: f32 },
}

impl FromStr for AssemblyArgument {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value == "asymmetric-unit" {
            return Ok(Self::AsymmetricUnit);
        }
        if let Some(id) = value.strip_prefix("biological:") {
            if id.is_empty() {
                return Err("biological assembly requires a non-empty ID".to_owned());
            }
            return Ok(Self::Biological(id.into()));
        }
        if let Some(text) = value.strip_prefix("crystal:") {
            let radius = positive_f32(text, "crystal radius")?;
            return Ok(Self::Crystal { radius });
        }
        Err("expected asymmetric-unit, biological:ID or crystal:RADIUS".to_owned())
    }
}

impl From<AssemblyArgument> for pdbiox::AssemblyChoice {
    fn from(value: AssemblyArgument) -> Self {
        match value {
            AssemblyArgument::AsymmetricUnit => Self::AsymmetricUnit,
            AssemblyArgument::Biological(id) => Self::Biological(id),
            AssemblyArgument::Crystal { radius } => Self::Crystal { radius },
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum ModelArgument {
    First,
    Index(u32),
    All,
    Ensemble,
}

impl FromStr for ModelArgument {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "first" => Ok(Self::First),
            "all" => Ok(Self::All),
            "ensemble" => Ok(Self::Ensemble),
            _ => {
                let Some(text) = value.strip_prefix("index:") else {
                    return Err("expected first, index:N, all or ensemble".to_owned());
                };
                text.parse::<u32>()
                    .map(Self::Index)
                    .map_err(|_| "model index must be an unsigned integer".to_owned())
            }
        }
    }
}

impl From<ModelArgument> for pdbiox::ModelChoice {
    fn from(value: ModelArgument) -> Self {
        match value {
            ModelArgument::First => Self::First,
            ModelArgument::Index(index) => Self::Index(index),
            ModelArgument::All => Self::All,
            ModelArgument::Ensemble => Self::Ensemble,
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) enum AltlocArgument {
    KeepAll,
    ConformerConsistent,
    Label(Box<str>),
    First,
    HighestOccupancyPerResidue,
    HighestOccupancyPerAtom,
}

impl FromStr for AltlocArgument {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "keep-all" => Ok(Self::KeepAll),
            "conformer-consistent" => Ok(Self::ConformerConsistent),
            "first" => Ok(Self::First),
            "highest-occupancy-per-residue" => Ok(Self::HighestOccupancyPerResidue),
            "highest-occupancy-per-atom" => Ok(Self::HighestOccupancyPerAtom),
            _ => {
                let Some(label) = value.strip_prefix("label:") else {
                    return Err("expected keep-all, conformer-consistent, label:ID, first, highest-occupancy-per-residue or highest-occupancy-per-atom".to_owned());
                };
                if label.is_empty() {
                    return Err("alternate-location label cannot be empty".to_owned());
                }
                Ok(Self::Label(label.into()))
            }
        }
    }
}

impl From<AltlocArgument> for pdbiox::AltlocPolicy {
    fn from(value: AltlocArgument) -> Self {
        match value {
            AltlocArgument::KeepAll => Self::KeepAll,
            AltlocArgument::ConformerConsistent => Self::ConformerConsistent,
            AltlocArgument::Label(label) => Self::Label(label),
            AltlocArgument::First => Self::First,
            AltlocArgument::HighestOccupancyPerResidue => Self::HighestOccupancyPerResidue,
            AltlocArgument::HighestOccupancyPerAtom => Self::HighestOccupancyPerAtom,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum NamespaceArgument {
    Label,
    Auth,
    Explicit,
}

impl From<NamespaceArgument> for pdbiox::Namespace {
    fn from(value: NamespaceArgument) -> Self {
        match value {
            NamespaceArgument::Label => Self::Label,
            NamespaceArgument::Auth => Self::Auth,
            NamespaceArgument::Explicit => Self::Explicit,
        }
    }
}

fn positive_f32(value: &str, label: &str) -> Result<f32, String> {
    match value.parse::<f32>() {
        Ok(number) if number.is_finite() && number > 0.0 => Ok(number),
        _ => Err(format!("{label} must be a positive finite number")),
    }
}
