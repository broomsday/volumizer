use numpy::{
    ndarray::{Array1, Array2},
    IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::{HashSet, VecDeque};
use std::i32;

const NATIVE_CONTRACT_VERSION: u32 = 1;
const OCCLUDED_DIMENSION_LIMIT: i32 = 4;

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
        return Err(PyValueError::new_err(
            "radius and coordinates must be finite",
        ));
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

fn validate_float_coordinate_array_shape(
    array: &PyReadonlyArray2<'_, f32>,
    name: &str,
) -> PyResult<()> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(PyValueError::new_err(format!(
            "{name} must have shape (N, 3), got {:?}",
            shape
        )));
    }
    Ok(())
}

#[pyfunction]
fn fibonacci_sphere_points_batch<'py>(
    py: Python<'py>,
    radius: f32,
    centers: PyReadonlyArray2<'_, f32>,
    samples: usize,
) -> PyResult<Bound<'py, PyArray2<f32>>> {
    if samples == 0 {
        return Err(PyValueError::new_err("samples must be > 0"));
    }
    if !radius.is_finite() {
        return Err(PyValueError::new_err("radius must be finite"));
    }
    validate_float_coordinate_array_shape(&centers, "centers")?;

    let centers_view = centers.as_array();
    let num_centers = centers_view.shape()[0];
    if num_centers == 0 {
        let array = Array2::from_shape_vec((0, 3), Vec::<f32>::new())
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        return Ok(array.into_pyarray_bound(py));
    }

    // Precompute sphere offsets once, then translate for each center.
    let golden_ratio = (1.0_f32 + 5.0_f32.sqrt()) / 4.0_f32;
    let mut offsets = Vec::with_capacity(samples * 3);
    for i in 0..samples {
        let i_f = i as f32;
        let samples_f = samples as f32;

        let phi = (1.0_f32 - 2.0_f32 * (i_f + 0.5_f32) / samples_f).acos();
        let theta = (std::f32::consts::PI * i_f) / golden_ratio;

        offsets.push((theta.cos() * phi.sin()) * radius);
        offsets.push((theta.sin() * phi.sin()) * radius);
        offsets.push(phi.cos() * radius);
    }

    let total_points = num_centers
        .checked_mul(samples)
        .ok_or_else(|| PyValueError::new_err("centers * samples is too large"))?;
    let mut values = Vec::with_capacity(total_points * 3);
    for center_index in 0..num_centers {
        let center = [
            centers_view[[center_index, 0]],
            centers_view[[center_index, 1]],
            centers_view[[center_index, 2]],
        ];
        if !center[0].is_finite() || !center[1].is_finite() || !center[2].is_finite() {
            return Err(PyValueError::new_err(format!(
                "centers[{}] must be finite, got {:?}",
                center_index, center
            )));
        }

        for offset_index in 0..samples {
            let base = offset_index * 3;
            values.push(center[0] + offsets[base]);
            values.push(center[1] + offsets[base + 1]);
            values.push(center[2] + offsets[base + 2]);
        }
    }

    let array = Array2::from_shape_vec((total_points, 3), values)
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

fn validate_grid_dimensions(grid_dimensions: &PyReadonlyArray1<'_, i32>) -> PyResult<[i32; 3]> {
    let dims = grid_dimensions.as_slice()?;
    if dims.len() != 3 {
        return Err(PyValueError::new_err(format!(
            "grid_dimensions must have shape (3,), got length {}",
            dims.len()
        )));
    }
    if dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0 {
        return Err(PyValueError::new_err(
            "grid_dimensions values must all be > 0",
        ));
    }

    Ok([dims[0], dims[1], dims[2]])
}

fn validate_voxel_in_grid(
    voxel: [i32; 3],
    grid_dimensions: [i32; 3],
    array_name: &str,
    index: usize,
) -> PyResult<()> {
    if voxel[0] < 0
        || voxel[1] < 0
        || voxel[2] < 0
        || voxel[0] >= grid_dimensions[0]
        || voxel[1] >= grid_dimensions[1]
        || voxel[2] >= grid_dimensions[2]
    {
        return Err(PyValueError::new_err(format!(
            "{array_name}[{index}] is out of bounds for grid_dimensions {:?}: {:?}",
            grid_dimensions, voxel
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

    let mut surface_indices: HashSet<i32> = direct_surface_indices
        .union(&neighbor_surface_indices)
        .copied()
        .collect();
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

    let mut reference_set: HashSet<[i32; 3]> =
        HashSet::with_capacity(reference.shape()[0].saturating_mul(2));
    for reference_index in 0..reference.shape()[0] {
        reference_set.insert([
            reference[[reference_index, 0]],
            reference[[reference_index, 1]],
            reference[[reference_index, 2]],
        ]);
    }

    let mut neighbor_indices: Vec<i32> = Vec::new();
    for query_index in 0..query.shape()[0] {
        let x = query[[query_index, 0]];
        let y = query[[query_index, 1]];
        let z = query[[query_index, 2]];

        // A query voxel is a 6-neighbor when any axis-adjacent coordinate exists in reference.
        if reference_set.contains(&[x - 1, y, z])
            || reference_set.contains(&[x + 1, y, z])
            || reference_set.contains(&[x, y - 1, z])
            || reference_set.contains(&[x, y + 1, z])
            || reference_set.contains(&[x, y, z - 1])
            || reference_set.contains(&[x, y, z + 1])
        {
            neighbor_indices.push(query_index as i32);
        }
    }

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
fn get_exposed_and_buried_voxel_indices(
    py: Python<'_>,
    solvent_voxels: PyReadonlyArray2<'_, i32>,
    protein_voxels: PyReadonlyArray2<'_, i32>,
    grid_dimensions: PyReadonlyArray1<'_, i32>,
) -> PyResult<PyObject> {
    validate_voxel_array_shape(&solvent_voxels, "solvent_voxels")?;
    validate_voxel_array_shape(&protein_voxels, "protein_voxels")?;
    let grid_dims = validate_grid_dimensions(&grid_dimensions)?;

    let grid_x = grid_dims[0] as usize;
    let grid_y = grid_dims[1] as usize;
    let grid_z = grid_dims[2] as usize;

    let z_array_len = grid_x
        .checked_mul(grid_y)
        .ok_or_else(|| PyValueError::new_err("grid_dimensions are too large"))?;
    let y_array_len = grid_x
        .checked_mul(grid_z)
        .ok_or_else(|| PyValueError::new_err("grid_dimensions are too large"))?;
    let x_array_len = grid_y
        .checked_mul(grid_z)
        .ok_or_else(|| PyValueError::new_err("grid_dimensions are too large"))?;

    // For each plane, track only min/max coordinate in the 3rd dimension.
    let mut z_min = vec![i32::MAX; z_array_len];
    let mut z_max = vec![i32::MIN; z_array_len];
    let mut y_min = vec![i32::MAX; y_array_len];
    let mut y_max = vec![i32::MIN; y_array_len];
    let mut x_min = vec![i32::MAX; x_array_len];
    let mut x_max = vec![i32::MIN; x_array_len];

    let protein_view = protein_voxels.as_array();
    for protein_index in 0..protein_view.shape()[0] {
        let voxel = [
            protein_view[[protein_index, 0]],
            protein_view[[protein_index, 1]],
            protein_view[[protein_index, 2]],
        ];
        validate_voxel_in_grid(voxel, grid_dims, "protein_voxels", protein_index)?;
        let x = voxel[0] as usize;
        let y = voxel[1] as usize;
        let z = voxel[2] as usize;

        let z_idx = x * grid_y + y;
        z_min[z_idx] = z_min[z_idx].min(voxel[2]);
        z_max[z_idx] = z_max[z_idx].max(voxel[2]);

        let y_idx = x * grid_z + z;
        y_min[y_idx] = y_min[y_idx].min(voxel[1]);
        y_max[y_idx] = y_max[y_idx].max(voxel[1]);

        let x_idx = y * grid_z + z;
        x_min[x_idx] = x_min[x_idx].min(voxel[0]);
        x_max[x_idx] = x_max[x_idx].max(voxel[0]);
    }

    let solvent_view = solvent_voxels.as_array();
    let mut exposed_indices: Vec<i32> = Vec::new();
    let mut buried_indices: Vec<i32> = Vec::new();

    for solvent_index in 0..solvent_view.shape()[0] {
        let voxel = [
            solvent_view[[solvent_index, 0]],
            solvent_view[[solvent_index, 1]],
            solvent_view[[solvent_index, 2]],
        ];
        validate_voxel_in_grid(voxel, grid_dims, "solvent_voxels", solvent_index)?;

        let x = voxel[0] as usize;
        let y = voxel[1] as usize;
        let z = voxel[2] as usize;
        let mut occluded_dimensions = [0_i32; 6];

        let z_idx = x * grid_y + y;
        if z_min[z_idx] != i32::MAX {
            if z_min[z_idx] < voxel[2] {
                occluded_dimensions[4] = 1;
            }
            if z_max[z_idx] > voxel[2] {
                occluded_dimensions[5] = 1;
            }
        }

        let y_idx = x * grid_z + z;
        if y_min[y_idx] != i32::MAX {
            if y_min[y_idx] < voxel[1] {
                occluded_dimensions[2] = 1;
            }
            if y_max[y_idx] > voxel[1] {
                occluded_dimensions[3] = 1;
            }
        }

        let x_idx = y * grid_z + z;
        if x_min[x_idx] != i32::MAX {
            if x_min[x_idx] < voxel[0] {
                occluded_dimensions[0] = 1;
            }
            if x_max[x_idx] > voxel[0] {
                occluded_dimensions[1] = 1;
            }
        }

        let num_occluded: i32 = occluded_dimensions.iter().sum();
        let buried = num_occluded > OCCLUDED_DIMENSION_LIMIT
            || (num_occluded == OCCLUDED_DIMENSION_LIMIT
                && ((occluded_dimensions[0] == 0 && occluded_dimensions[1] == 0)
                    || (occluded_dimensions[2] == 0 && occluded_dimensions[3] == 0)
                    || (occluded_dimensions[4] == 0 && occluded_dimensions[5] == 0)));

        if buried {
            buried_indices.push(solvent_index as i32);
        } else {
            exposed_indices.push(solvent_index as i32);
        }
    }

    let result = PyDict::new_bound(py);
    result.set_item(
        "exposed_indices",
        Array1::from(exposed_indices).into_pyarray_bound(py),
    )?;
    result.set_item(
        "buried_indices",
        Array1::from(buried_indices).into_pyarray_bound(py),
    )?;
    Ok(result.into())
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
    let grid_dims = validate_grid_dimensions(&grid_dimensions)?;
    if min_num_voxels < 0 {
        return Err(PyValueError::new_err("min_num_voxels must be >= 0"));
    }

    let buried_view = buried_voxels.as_array();
    let exposed_view = exposed_voxels.as_array();

    let buried: Vec<[i32; 3]> = (0..buried_view.shape()[0])
        .map(|i| {
            [
                buried_view[[i, 0]],
                buried_view[[i, 1]],
                buried_view[[i, 2]],
            ]
        })
        .collect();
    let exposed: Vec<[i32; 3]> = (0..exposed_view.shape()[0])
        .map(|i| {
            [
                exposed_view[[i, 0]],
                exposed_view[[i, 1]],
                exposed_view[[i, 2]],
            ]
        })
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
            get_agglomerated_type(&component_index_set, &buried, &exposed, grid_dims_array)
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
    module.add_function(wrap_pyfunction!(fibonacci_sphere_points_batch, module)?)?;
    module.add_function(wrap_pyfunction!(get_neighbor_voxel_indices, module)?)?;
    module.add_function(wrap_pyfunction!(bfs_component_indices, module)?)?;
    module.add_function(wrap_pyfunction!(
        get_exposed_and_buried_voxel_indices,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(classify_buried_components, module)?)?;
    Ok(())
}
