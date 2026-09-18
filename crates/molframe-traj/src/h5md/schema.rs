//! H5MD schema vocabulary.

pub(crate) const H5MD_GROUP: &str = "h5md";
pub(crate) const CREATOR_GROUP: &str = "h5md/creator";
pub(crate) const VERSION_ATTRIBUTE: &str = "version";
pub(crate) const CREATOR_NAME_ATTRIBUTE: &str = "name";
pub(crate) const CREATOR_VERSION_ATTRIBUTE: &str = "version";
pub(crate) const H5MD_VERSION: [i32; 2] = [1, 1];

pub(crate) const POSITION: &str = "position";
pub(crate) const VELOCITY: &str = "velocity";
pub(crate) const FORCE: &str = "force";
pub(crate) const VALUE: &str = "value";
pub(crate) const STEP: &str = "step";
pub(crate) const TIME: &str = "time";
pub(crate) const BOX: &str = "box";
pub(crate) const EDGES: &str = "edges";
pub(crate) const UNIT_ATTRIBUTE: &str = "unit";
pub(crate) const DIMENSION_ATTRIBUTE: &str = "dimension";
pub(crate) const BOUNDARY_ATTRIBUTE: &str = "boundary";
pub(crate) const PERIODIC: &str = "periodic";

pub(crate) fn element_path(group: &str, element: &str, member: &str) -> String {
    format!("particles/{group}/{element}/{member}")
}

pub(crate) fn box_path(group: &str, member: &str) -> String {
    format!("{}/{member}", box_edges_path(group))
}

pub(crate) fn box_edges_path(group: &str) -> String {
    format!("particles/{group}/{BOX}/{EDGES}")
}

pub(crate) fn box_group_path(group: &str) -> String {
    format!("particles/{group}/{BOX}")
}
