//! CPU `DLPack` tensors with mutable-safe, independently owned storage.

use crate::ExportCost;
use pdbiox_core::Structure;
use std::ffi::c_void;
use std::ptr::NonNull;

const CPU_DEVICE: i32 = 1;
const FLOAT_DATA_TYPE: u8 = 2;
const COORDINATE_DIMENSIONS: i32 = 2;
const COORDINATE_WIDTH: i64 = 3;
const FLOAT32_BITS: u8 = 32;

/// `DLPack` device descriptor.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DLDevice {
    /// Device kind; coordinate snapshots use the `DLPack` CPU code.
    pub device_type: i32,
    /// Device ordinal; CPU snapshots use device zero.
    pub device_id: i32,
}

/// `DLPack` scalar type descriptor.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DLDataType {
    /// Scalar kind; coordinate snapshots use the `DLPack` float code.
    pub code: u8,
    /// Bits per scalar lane.
    pub bits: u8,
    /// Vector lanes per element.
    pub lanes: u16,
}

/// ABI-compatible `DLPack` tensor view.
#[repr(C)]
#[derive(Debug)]
pub struct DLTensor {
    /// First byte of tensor storage.
    pub data: *mut c_void,
    /// Device carrying the storage.
    pub device: DLDevice,
    /// Number of dimensions.
    pub ndim: i32,
    /// Scalar representation.
    pub dtype: DLDataType,
    /// Mutable ABI pointer to dimension lengths owned by the manager context.
    pub shape: *mut i64,
    /// Null for compact row-major storage.
    pub strides: *mut i64,
    /// Byte displacement from `data`.
    pub byte_offset: u64,
}

/// ABI-compatible managed `DLPack` tensor.
#[repr(C)]
#[derive(Debug)]
pub struct DLManagedTensor {
    /// Tensor metadata and borrowed storage pointer.
    pub dl_tensor: DLTensor,
    /// Opaque owner for shape metadata and independent tensor storage.
    pub manager_ctx: *mut c_void,
    /// Releases tensor metadata and independent coordinate storage exactly once.
    pub deleter: Option<unsafe extern "C" fn(*mut DLManagedTensor)>,
}

struct ManagerContext {
    coordinates: Box<[f32]>,
    shape: Box<[i64; 2]>,
}

/// Owned producer-side handle for one `DLPack` managed tensor.
///
/// Coordinate export materialises one contiguous foreign-owned buffer because
/// `DLPack` has no read-only flag and consumers may legally mutate tensor data.
/// Dropping this handle releases that buffer unless ownership was transferred
/// with [`Self::into_raw`].
#[derive(Debug)]
pub struct DlpackTensor {
    managed: Option<NonNull<DLManagedTensor>>,
}

impl DlpackTensor {
    /// Exports model-zero coordinates as a compact `(atoms, 3)` CPU `float32`
    /// tensor using one contiguous copy isolated from the immutable snapshot.
    ///
    /// # Errors
    ///
    /// Returns an error when the topology does not describe one dense model.
    pub fn coordinates(structure: &Structure) -> Result<Self, DlpackError> {
        if structure.positions().len() != structure.atom_count() as usize {
            return Err(DlpackError::RaggedCoordinates);
        }
        let atoms = i64::from(structure.atom_count());
        let coordinates = structure
            .positions()
            .iter()
            .flat_map(|position| position.iter().copied())
            .collect::<Vec<_>>()
            .into_boxed_slice();
        let mut context = Box::new(ManagerContext {
            coordinates,
            shape: Box::new([atoms, COORDINATE_WIDTH]),
        });
        let data = context.coordinates.as_mut_ptr().cast::<c_void>();
        let shape = context.shape.as_mut_ptr();
        let manager_ctx = Box::into_raw(context).cast::<c_void>();
        let managed = Box::new(DLManagedTensor {
            dl_tensor: DLTensor {
                data,
                device: DLDevice {
                    device_type: CPU_DEVICE,
                    device_id: 0,
                },
                ndim: COORDINATE_DIMENSIONS,
                dtype: DLDataType {
                    code: FLOAT_DATA_TYPE,
                    bits: FLOAT32_BITS,
                    lanes: 1,
                },
                shape,
                strides: std::ptr::null_mut(),
                byte_offset: 0,
            },
            manager_ctx,
            deleter: Some(delete_managed),
        });
        Ok(Self {
            managed: NonNull::new(Box::into_raw(managed)),
        })
    }

    /// Cost classification for coordinate tensor export.
    ///
    /// This is [`ExportCost::Copy`] because `DLPack` cannot express read-only
    /// backing storage and pdbiox snapshots must remain immutable.
    #[must_use]
    pub const fn cost(&self) -> ExportCost {
        ExportCost::Copy
    }

    /// Borrows the managed tensor while producer ownership remains here.
    #[must_use]
    pub fn as_managed(&self) -> Option<&DLManagedTensor> {
        let managed = self.managed?;
        // SAFETY: `managed` originates from Box and remains owned by `self`.
        Some(unsafe { managed.as_ref() })
    }

    /// Transfers the managed tensor to a `DLPack` consumer.
    ///
    /// The consumer must call the published deleter exactly once. The producer
    /// no longer releases the tensor after this call.
    #[must_use]
    pub fn into_raw(mut self) -> *mut DLManagedTensor {
        match self.managed.take() {
            Some(managed) => managed.as_ptr(),
            None => std::ptr::null_mut(),
        }
    }
}

impl Drop for DlpackTensor {
    fn drop(&mut self) {
        let Some(managed) = self.managed.take() else {
            return;
        };
        // SAFETY: this handle retained producer ownership and calls once.
        unsafe { delete_managed(managed.as_ptr()) };
    }
}

/// Releases one producer allocation after `DLPack` ownership transfer.
///
/// # Safety
///
/// `managed` must be null or the pointer returned exactly once by
/// [`DlpackTensor::into_raw`]. `DLPack` consumers own that obligation.
unsafe extern "C" fn delete_managed(managed: *mut DLManagedTensor) {
    let Some(managed) = NonNull::new(managed) else {
        return;
    };
    // SAFETY: DLPack transfers this Box allocation to exactly one deleter call.
    let tensor = unsafe { Box::from_raw(managed.as_ptr()) };
    let context = tensor.manager_ctx.cast::<ManagerContext>();
    if let Some(context) = NonNull::new(context) {
        // SAFETY: manager_ctx was created by Box::into_raw beside this tensor.
        drop(unsafe { Box::from_raw(context.as_ptr()) });
    }
}

/// `DLPack` export failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum DlpackError {
    /// Topology rows do not form one rectangular coordinate block.
    #[error("DLPack coordinate export requires one dense model")]
    RaggedCoordinates,
}

#[cfg(test)]
#[path = "dlpack_tests.rs"]
mod tests;
