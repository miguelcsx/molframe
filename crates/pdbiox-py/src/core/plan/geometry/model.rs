//! Owned Python references for typed geometry operation nodes.

use pyo3::prelude::*;
use pyo3::types::PyAny;

#[derive(Debug)]
pub(crate) enum PyGeometryOperation {
    Centroid {
        positions: Py<PyAny>,
    },
    CentreOfMass {
        positions: Py<PyAny>,
        masses: Option<Py<PyAny>>,
    },
    RadiusOfGyration {
        positions: Py<PyAny>,
        masses: Option<Py<PyAny>>,
    },
    InertiaTensor {
        positions: Py<PyAny>,
        masses: Option<Py<PyAny>>,
    },
    PrincipalAxes {
        positions: Py<PyAny>,
        masses: Option<Py<PyAny>>,
        options: pdbiox::EigenOptions,
    },
    Asphericity {
        positions: Py<PyAny>,
        options: pdbiox::EigenOptions,
    },
    GyrationAxes {
        positions: Py<PyAny>,
        options: pdbiox::EigenOptions,
    },
    DistanceMatrix {
        positions: Py<PyAny>,
    },
    DistanceMatrixBetween {
        left: Py<PyAny>,
        right: Py<PyAny>,
    },
    Rmsf {
        frames: Py<PyAny>,
    },
}

impl PyGeometryOperation {
    pub(crate) fn clone_ref(&self, py: Python<'_>) -> Self {
        match self {
            Self::Centroid { positions } => Self::Centroid {
                positions: positions.clone_ref(py),
            },
            Self::CentreOfMass { positions, masses } => Self::CentreOfMass {
                positions: positions.clone_ref(py),
                masses: masses.as_ref().map(|value| value.clone_ref(py)),
            },
            Self::RadiusOfGyration { positions, masses } => Self::RadiusOfGyration {
                positions: positions.clone_ref(py),
                masses: masses.as_ref().map(|value| value.clone_ref(py)),
            },
            Self::InertiaTensor { positions, masses } => Self::InertiaTensor {
                positions: positions.clone_ref(py),
                masses: masses.as_ref().map(|value| value.clone_ref(py)),
            },
            Self::PrincipalAxes {
                positions,
                masses,
                options,
            } => Self::PrincipalAxes {
                positions: positions.clone_ref(py),
                masses: masses.as_ref().map(|value| value.clone_ref(py)),
                options: *options,
            },
            Self::Asphericity { positions, options } => Self::Asphericity {
                positions: positions.clone_ref(py),
                options: *options,
            },
            Self::GyrationAxes { positions, options } => Self::GyrationAxes {
                positions: positions.clone_ref(py),
                options: *options,
            },
            Self::DistanceMatrix { positions } => Self::DistanceMatrix {
                positions: positions.clone_ref(py),
            },
            Self::DistanceMatrixBetween { left, right } => Self::DistanceMatrixBetween {
                left: left.clone_ref(py),
                right: right.clone_ref(py),
            },
            Self::Rmsf { frames } => Self::Rmsf {
                frames: frames.clone_ref(py),
            },
        }
    }
}
