use super::parent_residue;

#[test]
fn every_variant_maps_to_its_parent_and_standard_names_are_unchanged() {
    for (variant, parent) in [
        ("HID", "HIS"),
        ("HIE", "HIS"),
        ("HIP", "HIS"),
        ("HSD", "HIS"),
        ("HSE", "HIS"),
        ("HSP", "HIS"),
        ("CYX", "CYS"),
        ("CYM", "CYS"),
        ("ASH", "ASP"),
        ("GLH", "GLU"),
        ("LYN", "LYS"),
        ("ALA", "ALA"),
        ("XYZ", "XYZ"),
    ] {
        assert_eq!(parent_residue(variant), parent);
    }
}
