use numpy::{
    ndarray::{Array1, Array2},
    IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::{HashSet, VecDeque};

const NATIVE_CONTRACT_VERSION: u32 = 1;

#[pyfunction]
fn contract_version() -> u32 {
    NATIVE_CONTRACT_VERSION
}

#[pyfunction]
fn backend_info() -> &'static str {
    "volumizer_native phase1 scaffold"
}

#[pyfunction]
fn fibonacci_sphere_points<'py>(
    py: Python<'py>,
    radius: f32,
    x: f32,
    y: f32,
    z: f32,
    samples: usize,
) -> PyResult<Bound<'py, PyArray2<f32>>> {
    if samples == 0 {
        return Err(PyValueError::new_err("samples must be > 0"));
    }
    if !radius.is_finite() || !x.is_finite() || !y.is_finite() || !z.is_finite() {
        return Err(PyValueError::new_err("radius and coordinates must be finite"));
    }

    // Keep this formula aligned with the existing Python implementation.
    let golden_ratio = (1.0_f32 + 5.0_f32.sqrt()) / 4.0_f32;
    let mut values = Vec::with_capacity(samples * 3);
    for i in 0..samples {
        let i_f = i as f32;
        let samples_f = samples as f32;

        let phi = (1.0_f32 - 2.0_f32 * (i_f + 0.5_f32) / samples_f).acos();
        let theta = (std::f32::consts::PI * i_f) / golden_ratio;

        values.push(x + (theta.cos() * phi.sin()) * radius);
        values.push(y + (theta.sin() * phi.sin()) * radius);
        values.push(z + phi.cos() * radius);
    }

    let array = Array2::from_shape_vec((samples, 3), values)
        .map_err(|err| PyValueError::new_err(err.to_string()))?;
    Ok(array.into_pyarray_bound(py))
}

fn is_neighbor_voxel(a: [i32; 3], b: [i32; 3]) -> bool {
    let mut sum_differences = 0_i32;
    for dimension in 0..3 {
        let difference = (a[dimension] - b[dimension]).abs();
        if difference > 1 {
            return false;
        }
        sum_differences += difference;
    }
    sum_differences == 1
}

fn validate_voxel_array_shape(array: &PyReadonlyArray2<'_, i32>, name: &str) -> PyResult<()> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err(format!(
            "{name} must have shape (N, 3), got {:?}",
            shape
        )));
    }
    Ok(())
}

fn bfs_component_from_set(voxels: &[[i32; 3]], searchable_indices: &HashSet<i32>) -> Vec<i32> {
    if searchable_indices.is_empty() {
        return Vec::new();
    }

    let start_index = *searchable_indices.iter().min().unwrap();
    let mut remaining_indices = searchable_indices.clone();
    remaining_indices.remove(&start_index);

    let mut queue: VecDeque<i32> = VecDeque::new();
    queue.push_back(start_index);

    let mut component_indices: Vec<i32> = vec![start_index];

    while let Some(current_index) = queue.pop_front() {
        let current_voxel = voxels[current_index as usize];

        let mut newly_found: Vec<i32> = Vec::new();
        for searched_index in remaining_indices.iter().copied() {
            let searched_voxel = voxels[searched_index as usize];
            if is_neighbor_voxel(current_voxel, searched_voxel) {
                newly_found.push(searched_index);
            }
        }

        for found_index in newly_found {
            if remaining_indices.remove(&found_index) {
                queue.push_back(found_index);
                component_indices.push(found_index);
            }
        }
    }

    component_indices.sort_unstable();
    component_indices
}

fn is_edge_voxel(voxel: [i32; 3], grid_dimensions: [i32; 3]) -> bool {
    voxel[0] == 0
        || voxel[1] == 0
        || voxel[2] == 0
        || voxel[0] == grid_dimensions[0] - 1
        || voxel[1] == grid_dimensions[1] - 1
        || voxel[2] == grid_dimensions[2] - 1
}

fn get_agglomerated_type(
    query_indices: &HashSet<i32>,
    buried_voxels: &[[i32; 3]],
    exposed_voxels: &[[i32; 3]],
    grid_dimensions: [i32; 3],
) -> (Vec<i32>, i8) {
    let mut direct_surface_indices: HashSet<i32> = HashSet::new();
    for query_index in query_indices.iter().copied() {
        let query_voxel = buried_voxels[query_index as usize];
        if is_edge_voxel(query_voxel, grid_dimensions) {
            direct_surface_indices.insert(query_index);
            continue;
        }

        for exposed_voxel in exposed_voxels.iter().copied() {
            if is_neighbor_voxel(query_voxel, exposed_voxel) {
                direct_surface_indices.insert(query_index);
                break;
            }
        }
    }

    // cavity
    if direct_surface_indices.is_empty() {
        return (Vec::new(), 1);
    }

    let mut neighbor_surface_indices: HashSet<i32> = HashSet::new();
    for query_index in query_indices
        .iter()
        .copied()
        .filter(|idx| !direct_surface_indices.contains(idx))
    {
        let query_voxel = buried_voxels[query_index as usize];
        for surface_index in direct_surface_indices.iter().copied() {
            let surface_voxel = buried_voxels[surface_index as usize];
            if is_neighbor_voxel(surface_voxel, query_voxel) {
                neighbor_surface_indices.insert(query_index);
                break;
            }
        }
    }

    let mut surface_indices: HashSet<i32> = direct_surface_indices.union(&neighbor_surface_indices).copied().collect();
    let single_surface_indices = bfs_component_from_set(buried_voxels, &surface_indices);
    if single_surface_indices.len() < surface_indices.len() {
        for index in single_surface_indices {
            surface_indices.remove(&index);
        }
        let second_surface = bfs_component_from_set(buried_voxels, &surface_indices);
        if second_surface.len() < surface_indices.len() {
            let mut sorted_surface_indices: Vec<i32> = direct_surface_indices
                .union(&neighbor_surface_indices)
                .copied()
                .collect();
            sorted_surface_indices.sort_unstable();
            // hub
            return (sorted_surface_indices, 4);
        }

        let mut sorted_surface_indices: Vec<i32> = direct_surface_indices
            .union(&neighbor_surface_indices)
            .copied()
            .collect();
        sorted_surface_indices.sort_unstable();
        // pore
        return (sorted_surface_indices, 3);
    }

    let mut sorted_surface_indices: Vec<i32> = direct_surface_indices
        .union(&neighbor_surface_indices)
        .copied()
        .collect();
    sorted_surface_indices.sort_unstable();
    // pocket
    (sorted_surface_indices, 2)
}

#[pyfunction]
fn get_neighbor_voxel_indices<'py>(
    py: Python<'py>,
    query_voxels: PyReadonlyArray2<'_, i32>,
    reference_voxels: PyReadonlyArray2<'_, i32>,
) -> PyResult<Bound<'py, PyArray1<i32>>> {
    validate_voxel_array_shape(&query_voxels, "query_voxels")?;
    validate_voxel_array_shape(&reference_voxels, "reference_voxels")?;

    let query = query_voxels.as_array();
    let reference = reference_voxels.as_array();

    let mut neighbor_indices: Vec<i32> = Vec::new();
    for query_index in 0..query.shape()[0] {
        let query_voxel = [
            query[[query_index, 0]],
            query[[query_index, 1]],
            query[[query_index, 2]],
        ];
        for reference_index in 0..reference.shape()[0] {
            let reference_voxel = [
                reference[[reference_index, 0]],
                reference[[reference_index, 1]],
                reference[[reference_index, 2]],
            ];
            if is_neighbor_voxel(query_voxel, reference_voxel) {
                neighbor_indices.push(query_index as i32);
                break;
            }
        }
    }

    neighbor_indices.sort_unstable();
    Ok(Array1::from(neighbor_indices).into_pyarray_bound(py))
}

#[pyfunction]
fn bfs_component_indices<'py>(
    py: Python<'py>,
    voxels: PyReadonlyArray2<'_, i32>,
    searchable_indices: PyReadonlyArray1<'_, i32>,
) -> PyResult<Bound<'py, PyArray1<i32>>> {
    validate_voxel_array_shape(&voxels, "voxels")?;
    let voxels_view = voxels.as_array();

    let num_voxels = voxels_view.shape()[0] as i32;
    let mut remaining_indices: HashSet<i32> = HashSet::new();
    for index in searchable_indices.as_slice()?.iter().copied() {
        if index < 0 || index >= num_voxels {
            return Err(PyValueError::new_err(format!(
                "searchable index out of bounds: {} (num voxels: {})",
                index, num_voxels
            )));
        }
        remaining_indices.insert(index);
    }

    if remaining_indices.is_empty() {
        return Ok(Array1::from(Vec::<i32>::new()).into_pyarray_bound(py));
    }

    let start_index = *remaining_indices
        .iter()
        .min()
        .ok_or_else(|| PyValueError::new_err("failed to initialize BFS start index"))?;
    remaining_indices.remove(&start_index);

    let mut queue: VecDeque<i32> = VecDeque::new();
    queue.push_back(start_index);

    let mut component_indices: Vec<i32> = vec![start_index];

    while let Some(current_index) = queue.pop_front() {
        let current_voxel = [
            voxels_view[[current_index as usize, 0]],
            voxels_view[[current_index as usize, 1]],
            voxels_view[[current_index as usize, 2]],
        ];

        let mut newly_found: Vec<i32> = Vec::new();
        for searched_index in remaining_indices.iter().copied() {
            let searched_voxel = [
                voxels_view[[searched_index as usize, 0]],
                voxels_view[[searched_index as usize, 1]],
                voxels_view[[searched_index as usize, 2]],
            ];
            if is_neighbor_voxel(current_voxel, searched_voxel) {
                newly_found.push(searched_index);
            }
        }

        for found_index in newly_found {
            if remaining_indices.remove(&found_index) {
                queue.push_back(found_index);
                component_indices.push(found_index);
            }
        }
    }

    component_indices.sort_unstable();
    Ok(Array1::from(component_indices).into_pyarray_bound(py))
}

#[pyfunction]
fn classify_buried_components(
    py: Python<'_>,
    buried_voxels: PyReadonlyArray2<'_, i32>,
    exposed_voxels: PyReadonlyArray2<'_, i32>,
    grid_dimensions: PyReadonlyArray1<'_, i32>,
    min_num_voxels: i32,
    _voxel_size: f32,
) -> PyResult<PyObject> {
    validate_voxel_array_shape(&buried_voxels, "buried_voxels")?;
    validate_voxel_array_shape(&exposed_voxels, "exposed_voxels")?;
    let grid_dims = grid_dimensions.as_slice()?;
    if grid_dims.len() != 3 {
        return Err(PyValueError::new_err(format!(
            "grid_dimensions must have shape (3,), got length {}",
            grid_dims.len()
        )));
    }
    if min_num_voxels < 0 {
        return Err(PyValueError::new_err("min_num_voxels must be >= 0"));
    }

    let buried_view = buried_voxels.as_array();
    let exposed_view = exposed_voxels.as_array();

    let buried: Vec<[i32; 3]> = (0..buried_view.shape()[0])
        .map(|i| [buried_view[[i, 0]], buried_view[[i, 1]], buried_view[[i, 2]]])
        .collect();
    let exposed: Vec<[i32; 3]> = (0..exposed_view.shape()[0])
        .map(|i| [exposed_view[[i, 0]], exposed_view[[i, 1]], exposed_view[[i, 2]]])
        .collect();
    let grid_dims_array = [grid_dims[0], grid_dims[1], grid_dims[2]];

    let buried_indices: HashSet<i32> = (0..(buried.len() as i32)).collect();
    let mut agglomerated_indices: HashSet<i32> = HashSet::new();

    let mut component_type_codes: Vec<i8> = Vec::new();
    let mut component_offsets: Vec<i32> = vec![0];
    let mut surface_offsets: Vec<i32> = vec![0];
    let mut component_voxel_indices_flat: Vec<i32> = Vec::new();
    let mut component_surface_indices_flat: Vec<i32> = Vec::new();

    while agglomerated_indices.len() < buried_indices.len() {
        let remaining_indices: HashSet<i32> = buried_indices
            .iter()
            .copied()
            .filter(|idx| !agglomerated_indices.contains(idx))
            .collect();
        let component_indices = bfs_component_from_set(&buried, &remaining_indices);
        for index in component_indices.iter().copied() {
            agglomerated_indices.insert(index);
        }

        let component_index_set: HashSet<i32> = component_indices.iter().copied().collect();
        let (surface_indices, type_code) = if component_indices.len() <= min_num_voxels as usize {
            (Vec::new(), 0_i8) // occluded
        } else {
            get_agglomerated_type(
                &component_index_set,
                &buried,
                &exposed,
                grid_dims_array,
            )
        };

        component_type_codes.push(type_code);
        component_voxel_indices_flat.extend(component_indices.iter().copied());
        component_surface_indices_flat.extend(surface_indices.iter().copied());
        component_offsets.push(component_voxel_indices_flat.len() as i32);
        surface_offsets.push(component_surface_indices_flat.len() as i32);
    }

    let result = PyDict::new_bound(py);
    result.set_item(
        "component_type_codes",
        Array1::from(component_type_codes).into_pyarray_bound(py),
    )?;
    result.set_item(
        "component_offsets",
        Array1::from(component_offsets).into_pyarray_bound(py),
    )?;
    result.set_item(
        "surface_offsets",
        Array1::from(surface_offsets).into_pyarray_bound(py),
    )?;
    result.set_item(
        "component_voxel_indices_flat",
        Array1::from(component_voxel_indices_flat).into_pyarray_bound(py),
    )?;
    result.set_item(
        "component_surface_indices_flat",
        Array1::from(component_surface_indices_flat).into_pyarray_bound(py),
    )?;

    Ok(result.into())
}

#[pymodule]
fn volumizer_native(_py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(contract_version, module)?)?;
    module.add_function(wrap_pyfunction!(backend_info, module)?)?;
    module.add_function(wrap_pyfunction!(fibonacci_sphere_points, module)?)?;
    module.add_function(wrap_pyfunction!(get_neighbor_voxel_indices, module)?)?;
    module.add_function(wrap_pyfunction!(bfs_component_indices, module)?)?;
    module.add_function(wrap_pyfunction!(classify_buried_components, module)?)?;
    Ok(())
}
