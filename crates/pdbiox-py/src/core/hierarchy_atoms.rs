use super::PyResidueAtoms;
use crate::atom::PyAtom;
use crate::errors::{index_error, key_error};
use crate::index::normalise_index;
use pdbiox::AtomRef;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

#[pymethods]
impl PyResidueAtoms {
    fn __len__(&self) -> usize {
        self.inner
            .residue(self.residue)
            .map_or(0, |residue| residue.atoms().count())
    }

    fn __getitem__(&self, key: &Bound<'_, PyAny>) -> PyResult<PyAtom> {
        let residue_key =
            isize::try_from(self.residue.get()).map_err(|_| index_error(key.py(), isize::MAX))?;
        let residue = self
            .inner
            .residue(self.residue)
            .ok_or_else(|| index_error(key.py(), residue_key))?;
        let atom = if let Ok(index) = key.extract::<isize>() {
            let Some(position) = normalise_index(index, residue.atoms().count()) else {
                return Err(index_error(key.py(), index));
            };
            residue.atom_at(position)
        } else if let Ok(name) = key.extract::<String>() {
            residue.atom(&name)
        } else {
            return Err(PyTypeError::new_err(
                "atom key must be an integer or string",
            ));
        };
        atom.map_or_else(
            || Err(key_error(key.py(), &key.str()?.to_string())),
            |atom| Ok(PyAtom::new(self.inner.clone(), AtomRef::index(atom))),
        )
    }
}
