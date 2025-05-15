use log::debug;
use num_complex::Complex32;
use phastft::fft_32_with_opts_and_plan;
use phastft::options::Options;
use phastft::planner::{Direction, Planner32};
use pyo3::prelude::*;

pub fn bit_reverse(val: usize, n_bits: usize) -> usize {
    let mut reversed = 0;
    match n_bits {
        0 => panic!("Can't have 0 bits!"),
        1 => reversed = val,
        _ => {
            for i in 0..n_bits {
                reversed |= ((val & (1 << i)) >> i) << (n_bits - i - 1);
            }
        }
    };
    reversed
}

#[pyfunction]
pub fn coarse_freq_correction_psk(
    iq_data: Vec<Complex32>,
    sample_rate: f32,
    psk_order: u32  // 1 for BPSK, 2 for QPSK, etc.
) -> (Vec<Complex32>, f32, usize) {
    // Taken from https://pysdr.org/content/sync.html#coarse-frequency-synchronization
    let next_p2 = 2usize.pow((iq_data.len() as f32).log2().ceil() as u32);
    println!("Padding input data length from {:?} to next power of 2 {:?} for FFT", iq_data.len(), next_p2);
    let mut reals: Vec<f32> = Vec::with_capacity(iq_data.len());
    let mut imags: Vec<f32> = Vec::with_capacity(iq_data.len());
    let mut corr_iq_data: Vec<Complex32> = iq_data.clone();
    for i in 0..next_p2 {
        if i >= iq_data.len() {
            reals.push(0.0);
            imags.push(0.0);
        } else {
            let c = iq_data[i].powu(psk_order + 1);
            reals.push(c.re);
            imags.push(c.im);
        }
    }
    let mut planner = Planner32::new(reals.len(), Direction::Forward);
    let opts = Options::default();
    fft_32_with_opts_and_plan(&mut reals, &mut imags, &opts, &mut planner);
    let psd: Vec<Complex32> = reals.iter().zip(imags.iter()).map(|(&r, &i)| Complex32::new(r, i)).collect();
    // let n_bits = (iq_data.len() as f64).log2() as usize;
    // for i in 0..iq_data.len() / 2 {
    //     let swap_i = bit_reverse(i, n_bits);
    //     psd.swap(i, swap_i);
    // }
    let mut max_freq_energy: f32 = 0.0;
    let mut offset_index: usize = 0;
    for i in 0..psd.len() {
        if psd[i].norm() > max_freq_energy {
            offset_index = i;
            max_freq_energy = psd[i].norm();
        }
    }
    let shifted_index = match offset_index >= next_p2 / 2 {
        true => offset_index - (next_p2 / 2),
        false => offset_index + (next_p2 / 2),
    };
    let freq_spacing = sample_rate / (next_p2 as f32);
    let freq_offset = (-sample_rate / 2.0) + (shifted_index as f32 * freq_spacing) - 0.5;
    println!("Offset index {:?} - Max energy {:?} - Freq offset: {:?}", offset_index, max_freq_energy, freq_offset);
    let ts: f32 = 1.0 / sample_rate;
    for i in 0..iq_data.len() {
        let shift = Complex32::new(0.0, -1.0 * (psk_order as f32) * core::f32::consts::PI * freq_offset * (i as f32) * ts).exp();
        corr_iq_data[i] = iq_data[i] * shift;
    }
    (corr_iq_data, freq_offset, offset_index)
}

#[pyfunction]
pub fn phastfft(iq_data: Vec<Complex32>) -> Vec<Complex32> {
    let mut reals: Vec<f32> = Vec::with_capacity(iq_data.len());
    let mut imags: Vec<f32> = Vec::with_capacity(iq_data.len());
    for iq in iq_data.iter() {
        reals.push(iq.re);
        imags.push(iq.im);
    }
    let mut planner = Planner32::new(reals.len(), Direction::Forward);
    let opts = Options::default();
    fft_32_with_opts_and_plan(&mut reals, &mut imags, &opts, &mut planner);
    reals.iter()
        .zip(imags.iter())
        .map(|(&r, &i)| Complex32::new(r, i))
        .collect()
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
    m.add_function(wrap_pyfunction!(phastfft, m)?)?;
    m.add_function(wrap_pyfunction!(costas_loop, m)?)?;
    Ok(())
}
