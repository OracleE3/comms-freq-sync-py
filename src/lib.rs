use log::debug;
use num_complex::Complex32;
use phastft::fft_32_with_opts_and_plan;
use phastft::options::Options;
use phastft::planner::{Direction, Planner32};
use pyo3::prelude::*;

#[pyfunction]
pub fn coarse_freq_correction_psk(
    iq_data: Vec<Complex32>,
    sample_rate: f32,
    psk_order: u32
) -> (Vec<Complex32>, f32) {
    // Taken from https://pysdr.org/content/sync.html#coarse-frequency-synchronization
    let mut reals: Vec<f32> = Vec::with_capacity(iq_data.len());
    let mut imags: Vec<f32> = Vec::with_capacity(iq_data.len());
    let mut corr_iq_data: Vec<Complex32> = iq_data.clone();
    for iq in iq_data.iter() {
        let c = iq.powu(psk_order);
        reals.push(c.re);
        imags.push(c.im);
    }
    let mut planner = Planner32::new(reals.len(), Direction::Forward);
    let opts = Options::default();
    fft_32_with_opts_and_plan(&mut reals, &mut imags, &opts, &mut planner);
    let offset_index = reals
        .iter()
        .zip(imags.iter())
        .enumerate()
        .map(|(index, (&r, &i))| (index, Complex32::new(r, i).norm()))
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .expect("Unable to get max value of PSD (freq offset)");
    // let freq_offset = (-sample_rate / 2.0) + ((offset_index as f32) * sample_rate / iq_data.len() as f32);
    let freq_offset = (offset_index as f32) * sample_rate / iq_data.len() as f32;
    debug!("Offset index {:?} - Freq offset: {:?}", offset_index, freq_offset);
    let ts: f32 = 1.0 / sample_rate;
    for i in 0..iq_data.len() {
        let shift = Complex32::new(0.0, -2.0 * (psk_order as f32) * core::f32::consts::PI * freq_offset * (i as f32) * ts).exp();
        corr_iq_data[i] = iq_data[i] * shift;
    }
    (corr_iq_data, freq_offset)
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
#[pyo3(name = "comms_freq_sync_py")]
fn src(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(coarse_freq_correction_psk, m)?)?;
    m.add_function(wrap_pyfunction!(costas_loop, m)?)?;
    Ok(())
}
