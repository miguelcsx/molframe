use super::{FIELD_COUNT, FIELD_NAMES, field_index};
use molframe_cif::Field;

#[test]
fn every_column_slot_is_the_position_of_the_field_it_names() {
    assert_eq!(FIELD_NAMES.len(), FIELD_COUNT);
    for (slot, name) in FIELD_NAMES.iter().enumerate() {
        assert_eq!(field_index(name), Some(slot), "{name}");
        let field = Field::from_item(name).unwrap_or_else(|| panic!("{name} is a field"));
        assert_eq!(field.position(), slot, "{name}");
    }
}
