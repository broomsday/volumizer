use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

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
) -> PyResult<&'py PyArray2<f32>> {
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
    Ok(array.into_pyarray(py))
}

#[pymodule]
fn volumizer_native(_py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(contract_version, module)?)?;
    module.add_function(wrap_pyfunction!(backend_info, module)?)?;
    module.add_function(wrap_pyfunction!(fibonacci_sphere_points, module)?)?;
    Ok(())
}
