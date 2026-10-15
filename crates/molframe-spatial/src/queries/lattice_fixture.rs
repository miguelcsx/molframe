//! Dense cubic lattice shared by the query tests.

/// Unit-spaced `side`-cubed lattice, dense enough to select the cell-list backend.
pub(crate) fn lattice(side: u16) -> Vec<[f32; 3]> {
    let mut positions = Vec::new();
    for x in 0..side {
        for y in 0..side {
            for z in 0..side {
                positions.push([f32::from(x), f32::from(y), f32::from(z)]);
            }
        }
    }
    positions
}
