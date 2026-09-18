//! Typed selection intermediate representation.

use molframe_chem::SmartsPattern;

/// A selectable structure column.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[non_exhaustive]
pub(crate) enum Column {
    Index,
    ByNumber,
    AtomSiteId,
    ResidueIndex,
    ChainIndex,
    ModelIndex,
    Chain,
    LabelChain,
    AuthChain,
    ResidueId,
    LabelResidueId,
    AuthResidueId,
    ResidueName,
    LabelResidueName,
    AuthResidueName,
    AtomName,
    LabelAtomName,
    AuthAtomName,
    AlternateLocation,
    Model,
    Entity,
    EntityType,
    Element,
    SegmentId,
    InsertionCode,
    RecordType,
    X,
    Y,
    Z,
    Mass,
    Charge,
    FormalCharge,
    Radius,
    BFactor,
    Occupancy,
    Plddt,
    Pae,
}

impl Column {
    #[must_use]
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        Some(match name.to_ascii_lowercase().as_str() {
            "index" => Self::Index,
            "bynum" => Self::ByNumber,
            "id" => Self::AtomSiteId,
            "resindex" => Self::ResidueIndex,
            "chainindex" => Self::ChainIndex,
            "modelindex" => Self::ModelIndex,
            "chain" => Self::Chain,
            "label_chain" => Self::LabelChain,
            "auth_chain" => Self::AuthChain,
            "resid" => Self::ResidueId,
            "label_resid" => Self::LabelResidueId,
            "auth_resid" => Self::AuthResidueId,
            "resname" => Self::ResidueName,
            "label_resname" => Self::LabelResidueName,
            "auth_resname" => Self::AuthResidueName,
            "name" => Self::AtomName,
            "label_name" => Self::LabelAtomName,
            "auth_name" => Self::AuthAtomName,
            "altloc" => Self::AlternateLocation,
            "model" => Self::Model,
            "entity" => Self::Entity,
            "entity_type" => Self::EntityType,
            "element" => Self::Element,
            "segid" => Self::SegmentId,
            "icode" => Self::InsertionCode,
            "record_type" => Self::RecordType,
            "x" => Self::X,
            "y" => Self::Y,
            "z" => Self::Z,
            "mass" => Self::Mass,
            "charge" => Self::Charge,
            "formalcharge" => Self::FormalCharge,
            "radius" => Self::Radius,
            "bfactor" | "tempfactor" => Self::BFactor,
            "occupancy" => Self::Occupancy,
            "plddt" => Self::Plddt,
            "pae" => Self::Pae,
            _ => return None,
        })
    }

    #[must_use]
    pub(crate) const fn is_numeric(self) -> bool {
        matches!(
            self,
            Self::Index
                | Self::ByNumber
                | Self::AtomSiteId
                | Self::ResidueIndex
                | Self::ChainIndex
                | Self::ModelIndex
                | Self::Model
                | Self::X
                | Self::Y
                | Self::Z
                | Self::Mass
                | Self::Charge
                | Self::FormalCharge
                | Self::Radius
                | Self::BFactor
                | Self::Occupancy
                | Self::Plddt
                | Self::Pae
        )
    }
}

/// A numeric comparison operator.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) enum Operator {
    Less,
    LessEqual,
    Greater,
    GreaterEqual,
    Equal,
    NotEqual,
}

impl Operator {
    #[must_use]
    pub(crate) const fn flipped(self) -> Self {
        match self {
            Self::Less => Self::Greater,
            Self::LessEqual => Self::GreaterEqual,
            Self::Greater => Self::Less,
            Self::GreaterEqual => Self::LessEqual,
            Self::Equal => Self::Equal,
            Self::NotEqual => Self::NotEqual,
        }
    }
}

/// A built-in chemistry or topology shorthand.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[non_exhaustive]
pub(crate) enum Macro {
    Protein,
    Backbone,
    Sidechain,
    Nucleic,
    NucleicBackbone,
    NucleicBase,
    NucleicSugar,
    Water,
    Ion,
    Lipid,
    Saccharide,
    Hetero,
    Hydrogen,
    Heavy,
    Polymer,
    Ligand,
    Aromatic,
}

impl Macro {
    #[must_use]
    pub(crate) fn from_name(name: &str) -> Option<Self> {
        Some(match name.to_ascii_lowercase().as_str() {
            "protein" => Self::Protein,
            "backbone" => Self::Backbone,
            "sidechain" => Self::Sidechain,
            "nucleic" => Self::Nucleic,
            "nucleicbackbone" => Self::NucleicBackbone,
            "nucleicbase" => Self::NucleicBase,
            "nucleicsugar" => Self::NucleicSugar,
            "water" => Self::Water,
            "ion" => Self::Ion,
            "lipid" => Self::Lipid,
            "saccharide" => Self::Saccharide,
            "hetero" => Self::Hetero,
            "hydrogen" => Self::Hydrogen,
            "heavy" => Self::Heavy,
            "polymer" => Self::Polymer,
            "ligand" => Self::Ligand,
            "aromatic" => Self::Aromatic,
            _ => return None,
        })
    }
}

/// The relation expanded by `same ... as`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
#[non_exhaustive]
pub(crate) enum SameKey {
    Residue,
    Chain,
    Model,
    Entity,
    Fragment,
    Segment,
    Column(Column),
}

/// A spatial selection operation.
#[derive(Clone, PartialEq, Debug)]
#[non_exhaustive]
pub(crate) enum GeometricExpr {
    Within {
        radius: f32,
        target: Box<Expr>,
    },
    Beyond {
        radius: f32,
        target: Box<Expr>,
    },
    Around {
        radius: f32,
        target: Box<Expr>,
    },
    SphereZone {
        radius: f32,
        target: Box<Expr>,
    },
    SphereLayer {
        inner: f32,
        outer: f32,
        target: Box<Expr>,
    },
    IsoLayer {
        inner: f32,
        outer: f32,
        target: Box<Expr>,
    },
    CylinderZone {
        radius: f32,
        z_max: f32,
        z_min: f32,
        target: Box<Expr>,
    },
    CylinderLayer {
        inner: f32,
        outer: f32,
        z_max: f32,
        z_min: f32,
        target: Box<Expr>,
    },
    Point {
        point: [f32; 3],
        radius: f32,
    },
}

/// A typed, storage-independent selection expression.
#[derive(Clone, PartialEq, Debug)]
#[non_exhaustive]
pub(crate) enum Expr {
    All,
    None,
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
    Not(Box<Self>),
    Global(Box<Self>),
    ByResidue(Box<Self>),
    Same {
        key: SameKey,
        target: Box<Self>,
    },
    Bonded {
        depth: u32,
        target: Box<Self>,
    },
    Geometric(GeometricExpr),
    Comparison {
        column: Column,
        operator: Operator,
        value: f64,
        absolute: bool,
    },
    Membership {
        column: Column,
        values: Vec<Box<str>>,
    },
    Group(Box<str>),
    Atom {
        segment: Box<str>,
        residue: i32,
        name: Box<str>,
    },
    Macro(Macro),
    Chirality(Box<str>),
    Smarts(SmartsPattern),
}
