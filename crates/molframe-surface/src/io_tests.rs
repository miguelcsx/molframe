use super::*;
use crate::SurfaceFace;

#[test]
fn obj_preserves_deterministic_indices_and_refuses_replacement() {
    let directory = match tempfile::tempdir() {
        Ok(directory) => directory,
        Err(error) => panic!("temporary directory failed: {error}"),
    };
    let path = directory.path().join("surface.obj");
    let mesh = IndexedSurfaceMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![SurfaceFace([0, 1, 2])],
    );
    if let Err(error) = write_obj(&path, &mesh) {
        panic!("OBJ write failed: {error}")
    }
    let text = match std::fs::read_to_string(&path) {
        Ok(text) => text,
        Err(error) => panic!("OBJ read failed: {error}"),
    };
    assert!(text.contains("f 1//1 2//2 3//3"));
    assert!(write_obj(path, &mesh).is_err());
}
