//! Deterministic tetrahedral polygonisation of an SES distance field.

use super::{SolventExcludedSurface, SurfaceTriangle};
use crate::SasaError;
use crate::cavity::Grid;
use crate::numeric::f64_to_f32;

const TETRAHEDRA: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 3, 2, 7],
    [0, 2, 6, 7],
    [0, 6, 4, 7],
    [0, 4, 5, 7],
    [0, 5, 1, 7],
];

/// Stack-only partition of tetrahedron vertices around the isosurface level.
struct TetraPartition {
    inside: [usize; 4],
    inside_len: usize,
    outside: [usize; 4],
    outside_len: usize,
}

/// Counts the exact topological triangle capacity without constructing output.
///
/// Degenerate geometric triangles may later be discarded, so this is a strict
/// upper bound and guarantees that polygonisation never reallocates.
pub(super) fn triangle_capacity(
    grid: &Grid,
    distance: &[f32],
    level: f32,
) -> Result<usize, SasaError> {
    let mut count = 0usize;
    for z in 0..grid.dims[2] - 1 {
        for y in 0..grid.dims[1] - 1 {
            for x in 0..grid.dims[0] - 1 {
                let values = cube_values(grid, distance, x, y, z);
                if !crosses_level(values, level) {
                    continue;
                }
                for tetrahedron in TETRAHEDRA {
                    count = count
                        .checked_add(tetrahedron_triangle_count(values, level, tetrahedron))
                        .ok_or(SasaError::GridDimensionsOverflow)?;
                }
            }
        }
    }
    Ok(count)
}

/// Polygonises the distance field at `level`.
///
/// Uniform cubes are rejected before Cartesian vertices are constructed.
/// Triangle area is accumulated as triangles are emitted, avoiding a second
/// full traversal of the output mesh.
pub(super) fn polygonise(
    grid: &Grid,
    distance: &[f32],
    level: f32,
    triangle_capacity: usize,
) -> Result<SolventExcludedSurface, SasaError> {
    let mut triangles = crate::workspace::empty_with_capacity(triangle_capacity)?;
    let mut area = 0.0;

    for z in 0..grid.dims[2] - 1 {
        for y in 0..grid.dims[1] - 1 {
            for x in 0..grid.dims[0] - 1 {
                let values = cube_values(grid, distance, x, y, z);

                if !crosses_level(values, level) {
                    continue;
                }

                let vertices = cube_vertices(grid, x, y, z);

                for tetrahedron in TETRAHEDRA {
                    polygonise_tetrahedron(
                        &vertices,
                        values,
                        level,
                        tetrahedron,
                        &mut triangles,
                        &mut area,
                    );
                }
            }
        }
    }

    Ok(SolventExcludedSurface { triangles, area })
}

/// Returns whether one cube contains vertices on both sides of `level`.
///
/// Runtime and auxiliary space are `O(8)`.
fn crosses_level(values: [f32; 8], level: f32) -> bool {
    let mut inside = false;
    let mut outside = false;

    for value in values {
        if value > level {
            inside = true;
        } else {
            outside = true;
        }

        if inside && outside {
            return true;
        }
    }

    false
}

/// Returns the eight Cartesian cube-centre vertices.
///
/// Only one `Grid::centre` call is needed; the remaining vertices are derived
/// by adding the constant grid step.
fn cube_vertices(grid: &Grid, x: usize, y: usize, z: usize) -> [[f64; 3]; 8] {
    let [x0, y0, z0] = grid.centre(x, y, z);
    let x1 = x0 + grid.step;
    let y1 = y0 + grid.step;
    let z1 = z0 + grid.step;

    [
        [x0, y0, z0],
        [x1, y0, z0],
        [x0, y1, z0],
        [x1, y1, z0],
        [x0, y0, z1],
        [x1, y0, z1],
        [x0, y1, z1],
        [x1, y1, z1],
    ]
}

/// Returns scalar values at the eight cube vertices.
///
/// Missing internal distance entries become infinity rather than panicking.
fn cube_values(grid: &Grid, distance: &[f32], x: usize, y: usize, z: usize) -> [f32; 8] {
    let base = grid.index(x, y, z);
    let nx = grid.dims[0];
    let plane = nx * grid.dims[1];
    [
        distance_at(distance, base),
        distance_at(distance, base + 1),
        distance_at(distance, base + nx),
        distance_at(distance, base + nx + 1),
        distance_at(distance, base + plane),
        distance_at(distance, base + plane + 1),
        distance_at(distance, base + plane + nx),
        distance_at(distance, base + plane + nx + 1),
    ]
}

/// Reads one distance-field cell defensively.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn distance_at(distance: &[f32], index: usize) -> f32 {
    match distance.get(index).copied() {
        Some(value) => value,
        None => f32::INFINITY,
    }
}

/// Polygonises one tetrahedron without temporary heap vectors.
///
/// Vertex classification uses two fixed four-element arrays, eliminating the
/// two `Vec` allocations previously performed for every tetrahedron.
fn polygonise_tetrahedron(
    vertices: &[[f64; 3]; 8],
    values: [f32; 8],
    level: f32,
    tetrahedron: [usize; 4],
    output: &mut Vec<SurfaceTriangle>,
    total_area: &mut f64,
) {
    let partition = partition_tetrahedron(values, level, tetrahedron);

    match (partition.inside_len, partition.outside_len) {
        (1, 3) => emit_single_side(
            vertices,
            values,
            level,
            partition.inside[0],
            [
                partition.outside[0],
                partition.outside[1],
                partition.outside[2],
            ],
            output,
            total_area,
        ),
        (3, 1) => emit_single_side(
            vertices,
            values,
            level,
            partition.outside[0],
            [
                partition.inside[0],
                partition.inside[1],
                partition.inside[2],
            ],
            output,
            total_area,
        ),
        (2, 2) => emit_quad(
            vertices,
            values,
            level,
            [partition.inside[0], partition.inside[1]],
            [partition.outside[0], partition.outside[1]],
            output,
            total_area,
        ),
        _ => {}
    }
}

/// Returns how many triangles one tetrahedron split can emit.
fn tetrahedron_triangle_count(values: [f32; 8], level: f32, tetrahedron: [usize; 4]) -> usize {
    let partition = partition_tetrahedron(values, level, tetrahedron);
    match (partition.inside_len, partition.outside_len) {
        (1, 3) | (3, 1) => 1,
        (2, 2) => 2,
        _ => 0,
    }
}

/// Partitions tetrahedron vertices into `> level` and `<= level` sets.
///
/// Runtime and auxiliary space are `O(4)`.
fn partition_tetrahedron(values: [f32; 8], level: f32, tetrahedron: [usize; 4]) -> TetraPartition {
    let mut inside = [0usize; 4];
    let mut outside = [0usize; 4];
    let mut inside_len = 0usize;
    let mut outside_len = 0usize;

    for index in tetrahedron {
        if values[index] > level {
            inside[inside_len] = index;
            inside_len += 1;
        } else {
            outside[outside_len] = index;
            outside_len += 1;
        }
    }

    TetraPartition {
        inside,
        inside_len,
        outside,
        outside_len,
    }
}

/// Emits the single triangle arising from a one-versus-three tetrahedron split.
///
/// Runtime and auxiliary space are `O(1)`.
fn emit_single_side(
    vertices: &[[f64; 3]; 8],
    values: [f32; 8],
    level: f32,
    single: usize,
    others: [usize; 3],
    output: &mut Vec<SurfaceTriangle>,
    total_area: &mut f64,
) {
    push_triangle(
        output,
        total_area,
        [
            edge_intersection(vertices, values, single, others[0], level),
            edge_intersection(vertices, values, single, others[1], level),
            edge_intersection(vertices, values, single, others[2], level),
        ],
    );
}

/// Emits the two triangles arising from a two-versus-two tetrahedron split.
///
/// Runtime and auxiliary space are `O(1)`.
fn emit_quad(
    vertices: &[[f64; 3]; 8],
    values: [f32; 8],
    level: f32,
    inside: [usize; 2],
    outside: [usize; 2],
    output: &mut Vec<SurfaceTriangle>,
    total_area: &mut f64,
) {
    let ac = edge_intersection(vertices, values, inside[0], outside[0], level);
    let ad = edge_intersection(vertices, values, inside[0], outside[1], level);
    let bc = edge_intersection(vertices, values, inside[1], outside[0], level);
    let bd = edge_intersection(vertices, values, inside[1], outside[1], level);

    push_triangle(output, total_area, [ac, ad, bc]);
    push_triangle(output, total_area, [ad, bd, bc]);
}

/// Computes the isosurface intersection on one cube edge.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn edge_intersection(
    vertices: &[[f64; 3]; 8],
    values: [f32; 8],
    first: usize,
    second: usize,
    level: f32,
) -> [f64; 3] {
    interpolate(
        vertices[first],
        vertices[second],
        values[first],
        values[second],
        level,
    )
}

/// Linearly interpolates an isosurface point between two scalar samples.
///
/// Degenerate or non-finite scalar spans fall back to the edge midpoint.
fn interpolate(
    first: [f64; 3],
    second: [f64; 3],
    first_value: f32,
    second_value: f32,
    level: f32,
) -> [f64; 3] {
    let first_value = f64::from(first_value);
    let second_value = f64::from(second_value);
    let level = f64::from(level);
    let span = second_value - first_value;

    let fraction = if span.abs() <= f64::EPSILON || !span.is_finite() {
        0.5
    } else {
        ((level - first_value) / span).clamp(0.0, 1.0)
    };

    [
        first[0] + fraction * (second[0] - first[0]),
        first[1] + fraction * (second[1] - first[1]),
        first[2] + fraction * (second[2] - first[2]),
    ]
}

/// Emits one non-degenerate triangle and accumulates its area.
///
/// Degenerate triangles are discarded. Runtime and auxiliary space are `O(1)`.
fn push_triangle(output: &mut Vec<SurfaceTriangle>, total_area: &mut f64, vertices: [[f64; 3]; 3]) {
    let first = subtract(vertices[1], vertices[0]);
    let second = subtract(vertices[2], vertices[0]);
    let cross = cross(first, second);
    let length = vector_length(cross);

    if length <= f64::EPSILON {
        return;
    }

    let area = 0.5 * length;
    *total_area += area;

    output.push(SurfaceTriangle {
        vertices: vertices.map(|point| point.map(f64_to_f32)),
        normal: [
            f64_to_f32(cross[0] / length),
            f64_to_f32(cross[1] / length),
            f64_to_f32(cross[2] / length),
        ],
        area,
    });
}

/// Subtracts two three-dimensional vectors.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

/// Computes a three-dimensional cross product.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}

/// Returns the Euclidean length of a three-dimensional vector.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn vector_length(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}
