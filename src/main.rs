// -------------------------------------------------------------------------------------------------
//
//  bump - Computing the height of water for the case of a shallow water flow over a bump
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//!  bump - Computing the height of water for the case of a shallow water flow over a bump
//!
//!  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//!  All Rights Reserved.

#![allow(unused)]

//
// Application modules.
//
pub mod bump;
pub mod shallow_water;

//
// Import of internal modules
//
use crate::bump::Bump;

///
/// Main application driver.
///
fn main() {
    //
    // Print the application title on the screen.
    //
    println!(
        "\nBump - Computation of the height of water for the case\nof a shallow water flow over a bump\n\nVersion 3.0.0\n\n"
    );
    //
    // Create new Bump object.
    //
    let length: f64 = 25.0;
    let n_cells: usize = 200;
    let bump = Bump::new(length, n_cells);
    //
    // Subcritical flow over a bump.
    //
    let q_in: f64 = 4.42;
    let h_out: f64 = 2.0;
    let g: f64 = 9.81;
    let filename: &str = "subcritical_bump.dat";
    let screen_output: bool = true;
    bump.compute_subcritical_case(q_in, h_out, g, filename, screen_output);
    //
    // Transcritical no-shock flow over a bump.
    //
    let q_in: f64 = 1.53;
    let filename: &str = "transcritical_no_shock_bump.dat";
    bump.compute_transcritical_no_shock_case(q_in, g, filename, screen_output);
    //
    // Transcritical shock flow over a bump.
    //
    let q_in: f64 = 0.18;
    let h_out: f64 = 0.33;
    let filename: &str = "transcritical_shock_bump.dat";
    bump.compute_transcritical_shock_case(q_in, h_out, g, filename, screen_output);
}
