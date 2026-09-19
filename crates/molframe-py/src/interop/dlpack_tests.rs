use super::*;
use crate::structure::PyStructure;
use pyo3::types::PyCapsuleMethods;
use std::path::PathBuf;

#[test]
fn structure_publishes_a_standard_owned_dlpack_capsule() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    let structure = match molframe::read(path) {
        Ok(structure) => structure,
        Err(findings) => panic!("binding fixture failed to read: {findings:?}"),
    };
    let structure = PyStructure::new(structure);
    Python::initialize();
    Python::attach(|py| {
        assert_eq!(
            PyStructure::__dlpack_device__(),
            (CPU_DEVICE_TYPE, CPU_DEVICE_ID)
        );
        assert!(
            structure
                .__dlpack__(py, None, None, None, Some(false))
                .is_err()
        );
        let capsule = match structure.__dlpack__(py, None, None, None, None) {
            Ok(capsule) => capsule,
            Err(error) => panic!("DLPack capsule export failed: {error}"),
        };
        let pointer = match capsule.pointer_checked(Some(c"dltensor")) {
            Ok(pointer) => pointer.cast::<molframe::DLManagedTensor>(),
            Err(error) => panic!("DLPack capsule name or pointer is invalid: {error}"),
        };
        // SAFETY: the checked capsule pointer names a live DLManagedTensor.
        let tensor = unsafe { pointer.as_ref() };
        assert_eq!(tensor.dl_tensor.ndim, 2);
        assert_eq!(tensor.dl_tensor.device.device_type, CPU_DEVICE_TYPE);
        // Simulate the standard consumer ownership transition.
        let renamed =
            unsafe { pyo3::ffi::PyCapsule_SetName(capsule.as_ptr(), c"used_dltensor".as_ptr()) };
        assert_eq!(renamed, 0);
        if let Some(deleter) = tensor.deleter {
            // SAFETY: ownership was transferred away from the renamed capsule.
            unsafe { deleter(pointer.as_ptr()) };
        }
    });
}
