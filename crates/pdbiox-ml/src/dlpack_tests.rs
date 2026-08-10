use super::*;
use pdbiox_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "data_tensor\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 2 3\n\
ATOM 2 N N ALA A 1 4 5 6\n";

fn structure() -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn coordinates_publish_cpu_float32_shape_in_an_isolated_copy() {
    let structure = structure();
    let expected = structure.positions().as_ptr().cast::<f32>();
    let tensor = DlpackTensor::coordinates(&structure).expect("dense tensor");
    let managed = tensor.as_managed().expect("managed tensor");
    assert_ne!(managed.dl_tensor.data.cast_const().cast::<f32>(), expected);
    assert_eq!(managed.dl_tensor.ndim, 2);
    assert_eq!(managed.dl_tensor.dtype.code, 2);
    assert_eq!(managed.dl_tensor.dtype.bits, 32);
    assert!(managed.dl_tensor.strides.is_null());
    // SAFETY: shape has two entries retained by the live manager context.
    let shape = unsafe { std::slice::from_raw_parts(managed.dl_tensor.shape, 2) };
    assert_eq!(shape, &[2, 3]);
    assert_eq!(tensor.cost(), ExportCost::Copy);
    // SAFETY: the managed tensor owns six writable f32 values independently
    // from the immutable Structure snapshot.
    unsafe { managed.dl_tensor.data.cast::<f32>().write(9.0) };
    assert_eq!(structure.positions()[0][0].to_bits(), 1.0f32.to_bits());
}

#[test]
fn transferred_tensor_retains_its_copy_until_consumer_deletes_it() {
    let raw = {
        let structure = structure();
        DlpackTensor::coordinates(&structure)
            .expect("dense tensor")
            .into_raw()
    };
    assert!(!raw.is_null());
    // SAFETY: ownership was transferred above and remains live here.
    let managed = unsafe { &*raw };
    // SAFETY: two atoms times three coordinates are retained by manager_ctx.
    let coordinates =
        unsafe { std::slice::from_raw_parts(managed.dl_tensor.data.cast::<f32>(), 6) };
    assert_eq!(coordinates, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let deleter = managed.deleter.expect("consumer deleter");
    // SAFETY: the transferred consumer calls the deleter exactly once.
    unsafe { deleter(raw) };
}
