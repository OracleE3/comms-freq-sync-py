use std::f32::consts::PI;

use log::{debug};
use num::{complex::Complex32, Complex};
use phastft::{planner::Direction, fft_32};
use pyo3::pyfunction;

#[pyfunction]
pub fn coarse_freq_correction_psk(
    iq_data: Vec<Complex32>,
    sample_rate: f32,
    psk_order: usize
) -> Vec<Complex32> {
    // Taken from https://pysdr.org/content/sync.html#coarse-frequency-synchronization
    let mut reals: Vec<f32> = Vec::with_capacity(iq_data.len());
    let mut imags: Vec<f32> = Vec::with_capacity(iq_data.len());
    let mut corr_iq_data: Vec<Complex32> = Vec::with_capacity(iq_data.len());
    for iq in iq_data {
        let c = iq.pow(psk_order);
        reals.push(c.re());
        imags.push(c.im());
    }
    let psd = fft_32(&mut reals, &mut imags, Direction::Forward);
    let freq_offset = psd.max();
    debug!("Freq offset: {?}", freq_offset);
    let ts: f32 = 1.0 / sample_rate;
    for i in 0..corr_iq_data.len() {
        corr_iq_data[i] *= Complex32::new(0.0, 2.0 * psk_order * PI * ).exp()
    }
    corr_iq_data
}

#[pyfunction]
pub fn costas_loop(
    iq_data: Vec<Complex32>,
    sample_rate: f32,
    alpha: f32,
    beta: f32,
) -> (Vec<Complex32>, Vec<f32>) {
    // Costas loop
    let mut phase: f32 = 0.0;
    let mut freq: f32 = 0.0;
    let mut error: f32;

    let mut out: Vec<Complex32> = Vec::new();
    let mut freq_log: Vec<f32> = Vec::new();
    out.reserve(iq_data.len());
    freq_log.reserve(iq_data.len());

    for i in 0..iq_data.len() {
        // adjust the input sample by the inverse of the estimated phase offset
        out.push(iq_data[i] * Complex32::new(0.0, -phase).exp());
        // This is the error formula for 2nd order Costas Loop (e.g. for BPSK)
        error = out[i].re * out[i].im;

        // Advance the loop (recalc phase and freq offset)
        freq += beta * error;
        // convert from angular velocity to Hz for logging
        freq_log.push(freq * sample_rate / (2.0 * core::f32::consts::PI));
        phase += freq + (alpha * error);

        // Optional: Adjust phase so its always between 0 and 2pi, recall that phase wraps around every 2pi
        phase %= 2.0 * core::f32::consts::PI;
    }

    (out, freq_log)
}

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "_lib")]
fn src(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(costas_loop, m)?)?;
    Ok(())
}
