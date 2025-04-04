// -------------------------------------------------------------------------------------------------
//
//  bump - Computing the height of water for the case of a shallow water flow over a bump
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//
// Import of external crates
//
use num_complex::Complex64;

///
/// Compute the Bernoulli coefficients i.e. the coefficient of the cubic form
/// of the Bernoulli equation a*h^3 + b*h^2 + c*h + d = 0 for the shallow
/// water flow between an inflow and an outflow location.
///
/// - Arguments:
///   - 'q_in': inflow discharge (m²/s).
///   - 'h_out': outflow height (m).
///   - 'z_in': inflow elevation (m).
///   - 'z_out': outflow elevation (m).
///   - 'g': acceleration due to gravity (m/s²).
///
/// - Returns:
///   - The Bernoulli coefficients a, b, c, d.
///
pub fn bernoulli_coefficients(
    q_in: f64,
    h_out: f64,
    z_in: f64,
    z_out: f64,
    g: f64,
) -> (f64, f64, f64, f64) {
    let a = 1.0;
    let b = -(q_in * q_in / (2.0 * g * h_out * h_out) + h_out - (z_in - z_out));
    let c = 0.0;
    let d = q_in * q_in / (2.0 * g);

    (a, b, c, d)
}

///
/// Find the solution which makes sense for a particular case
/// given the three roots of the cubic equation.
///
/// - Arguments:
///   - 'roots': the roots of the cubic equation.
///   - 'h_near': The height of water at the previously computed location (m).
///
/// - Returns:
///   - The solution which makes sense.
///
pub fn find_solution(roots: &[Complex64], h_near: f64) -> f64 {
    //
    // Ignore the roots with imaginary parts and negative
    // real parts.
    //
    let mut real_roots: Vec<f64> = Vec::new();
    for x in roots.iter() {
        if x.im.abs() <= 1.0e-6 && x.re >= 0.0 {
            real_roots.push(x.re);
        }
    }
    //
    // Select the root which is closest to the previously
    // computed height of water.
    //
    let mut solution: f64 = 0.0;
    let mut min_diff: f64 = f64::MAX;
    for root in real_roots.iter() {
        let diff = (root - h_near).abs();
        if diff < min_diff {
            min_diff = diff;
            solution = *root;
        }
    }

    solution
}
