use serde::{Deserialize, Serialize};
use serde_bytes::ByteBuf;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct File {
    pub mmtf_version: String,
    pub mmtf_producer: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub structure_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit_cell: Option<Vec<f32>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub space_group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub experimental_methods: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolution: Option<f32>,
    pub num_bonds: i32,
    pub num_atoms: i32,
    pub num_groups: i32,
    pub num_chains: i32,
    pub num_models: i32,
    pub group_list: Vec<Group>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bond_atom_list: Option<ByteBuf>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bond_order_list: Option<ByteBuf>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bond_resonance_list: Option<ByteBuf>,
    pub x_coord_list: ByteBuf,
    pub y_coord_list: ByteBuf,
    pub z_coord_list: ByteBuf,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub b_factor_list: Option<ByteBuf>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub atom_id_list: Option<ByteBuf>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub alt_loc_list: Option<ByteBuf>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub occupancy_list: Option<ByteBuf>,
    pub group_id_list: ByteBuf,
    pub group_type_list: ByteBuf,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ins_code_list: Option<ByteBuf>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequence_index_list: Option<ByteBuf>,
    pub chain_id_list: ByteBuf,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chain_name_list: Option<ByteBuf>,
    pub groups_per_chain: Vec<i32>,
    pub chains_per_model: Vec<i32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub entity_list: Option<Vec<Entity>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct Group {
    pub formal_charge_list: Vec<i32>,
    pub atom_name_list: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub element_list: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bond_atom_list: Vec<i32>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bond_order_list: Vec<i32>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bond_resonance_list: Vec<i32>,
    #[serde(rename = "groupName")]
    pub name: String,
    pub single_letter_code: String,
    pub chem_comp_type: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct Entity {
    #[serde(default)]
    pub chain_index_list: Vec<i32>,
    #[serde(default)]
    pub description: String,
    #[serde(rename = "type", default)]
    pub kind: String,
    #[serde(default)]
    pub sequence: String,
}
