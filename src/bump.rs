// -------------------------------------------------------------------------------------------------
//
//  bump - Computing the height of water for the case of a shallow water flow over a bump
//
//  Copyright (c) 2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//
// Imports from the standard library.
//
use std::{io::Write, ops::Sub};

//
// Import of internal modules
//
use crate::shallow_water;

///
/// This struct contains the geometrical parameters of a Bump needed for the
/// numerical solution of the shallow water equations.
///
pub struct Bump {
    ///
    /// The length of the domain [m].
    ///
    length: f64,
    ///
    /// The number of cells to be used for the numerical discretisation.
    ///
    n_cells: usize,
    ///
    /// The x-coordinates of the cell centers [m].
    ///
    x: Vec<f64>,
    ///
    /// The height of the bump at the cell centers [m].
    ///
    z: Vec<f64>,
    ///
    /// The maximum height of the bump [m].
    ///
    z_max: f64,
    ///
    /// The location of the maximum height of the bump.
    ///
    n_max: usize,
}

//
// Implementation of the Bump struct.
//
impl Bump {
    ///
    /// Create a new Bump object and initialise its geometry.
    ///
    /// - Arguments:
    ///   - `length`: the length of the domain [m].
    ///   - `n_cells`: the number of cells to be used for the numerical discretisation.
    ///
    /// - Returns:
    ///   - A new Bump object.
    ///
    pub fn new(length: f64, n_cells: usize) -> Self {
        let mut bump = Self {
            length,
            n_cells,
            x: vec![0.0; n_cells + 1],
            z: vec![0.0; n_cells + 1],
            z_max: 0.0,
            n_max: 0,
        };
        //
        // Discretisation of the domain.
        //
        let dx = length / (n_cells as f64);
        for n in 0..=n_cells {
            bump.x[n] = (n as f64 - 0.5) * dx;
            bump.z[n] = Bump::bump_geometry(bump.x[n]);
        }
        //
        // Compute the maximum bump height and its location.
        //
        bump.z_max = 0.0;
        bump.n_max = 0;
        for n in 0..=n_cells {
            if bump.z[n] > bump.z_max {
                bump.z_max = bump.z[n];
                bump.n_max = n;
            }
        }
        bump.z_max = 0.2;

        bump
    }

    ///
    /// Compute the geometry of the bump.
    ///
    /// - Arguments:
    ///  - `x`: the x-coordinate [m].
    ///
    /// - Returns:
    ///   - The height of the bump for the given x-coordinate [m].
    ///
    #[inline]
    fn bump_geometry(x: f64) -> f64 {
        f64::max(0.0, 0.2 - 0.05 * (x - 10.0) * (x - 10.0))
    }

    ///
    /// Computes the Rankine-Hugoniot condition for shock conditions.
    ///
    /// - Arguments:
    ///   - `h_plus`: the height of the water on the right side of the shock [m].
    ///   - `h_minus`: the height of the water on the left side of the shock [m].
    ///   - `q`: flow rate of water [m²/s].
    ///   - `g`: the acceleration due to gravity [m/s²].
    ///
    /// - Returns:
    ///   - The value of the Rankine-Hugoniot condition which should be zero
    ///     if the values supplied as arguments satisfy it.
    ///
    #[inline]
    fn rankine_hugoniot(h_plus: f64, h_minus: f64, q: f64, g: f64) -> f64 {
        (q * q * (1.0 / h_plus - 1.0 / h_minus) + 0.5 * g * (h_plus - h_minus) * (h_plus + h_minus))
            .abs()
    }

    ///
    /// Computes the subcritical case and outputs the results.
    ///
    /// - Arguments:
    ///   - `q_in`: the inflow discharge [m²/s].
    ///   - `h_out`: the outflow height [m].
    ///   - `g`: the acceleration due to gravity [m/s²].
    ///   - `filename`: the name of the file to output the results.
    ///   - `screen_output`: whether to output the results to the screen.
    ///
    pub fn compute_subcritical_case(
        &self,
        q_in: f64,
        h_out: f64,
        g: f64,
        filename: &str,
        screen_output: bool,
    ) {
        //
        // Subcritical flow over a bump.
        //
        if screen_output {
            println!();
            println!("Subcritical Flow Over a Bump");
            println!("==============================================");
            println!("x (m)     Z (m)     height (m)");
        }
        //
        // Open a file to save the height of water.
        //
        let test_file = std::fs::File::create(filename);
        let mut subcritical_file = match test_file {
            Ok(file) => file,
            Err(_) => {
                println!("Bump - Error opening file {}", filename);
                return;
            }
        };
        writeln!(subcritical_file, "# x (m)    height (m)");
        //
        // Iteration from the end of the field to the beginning
        // with application of the Bernoulli equation at each cell.
        //
        let mut h_water = vec![0.0; self.n_cells + 1];
        h_water[self.n_cells] = h_out;

        for n in (0..self.n_cells).rev() {
            //
            // The coefficients of the Bernoulli cubic equation to solve
            // a*h^3 + b*h^2 + c*h + d = 0.
            //
            let coeffs = shallow_water::bernoulli_coefficients(q_in, h_out, self.z[n], 0.0, g);
            //
            // The bump height solutions.
            //
            let height = shallow_water::solve_cubic_equation(&coeffs);
            //
            // Find the solution which makes sense for this case.
            //
            h_water[n] = shallow_water::find_solution(&height, h_water[n + 1]);
        }
        //
        // Use zip to iterate over corresponding elements in both vectors
        // and update the water height.
        //
        h_water.iter_mut().zip(&self.z).for_each(|(h, &z)| {
            *h += z;
        });
        //
        // Output the results to the screen and a file.
        //
        for n in 1..=self.n_cells {
            let x = self.x[n];
            let z = self.z[n];
            let h = h_water[n];
            //
            // Output the results to the screen.
            //
            if screen_output {
                println!("{:.6}  {:.6}  {:.6}", x, z, h);
            }
            //
            // Output the results to a file.
            //
            writeln!(subcritical_file, "{:.6}   {:.6}", x, h).unwrap();
        }

        subcritical_file.flush().ok();
    }

    ///
    /// Computes the transcritical no-shock case and outputs the results.
    ///
    /// - Arguments:
    ///   - `q_in`: the inflow discharge [m²/s].
    ///   - `g`: the acceleration due to gravity [m/s²].
    ///   - `filename`: the name of the file to output the results.
    ///   - `screen_output`: whether to output the results to the screen.
    ///
    pub fn compute_transcritical_no_shock_case(
        &self,
        q_in: f64,
        g: f64,
        filename: &str,
        screen_output: bool,
    ) {
        //
        // Trancritical no-shock flow over a bump.
        //
        if screen_output {
            println!();
            println!("Transcritical no-Shock Flow Over a Bump");
            println!("==============================================");
            println!("x (m)     Z (m)     height (m)");
        }
        //
        // Open a file to save the height of water.
        //
        let test_file = std::fs::File::create(filename);
        let mut transcritical_no_shock_file = match test_file {
            Ok(file) => file,
            Err(_) => {
                println!("Bump - Error opening file {}", filename);
                return;
            }
        };
        writeln!(transcritical_no_shock_file, "# x (m)    height (m)");
        //
        // Discretisation.
        //
        let epsilon = 1.0 / self.n_cells as f64;
        //
        // Water height at the top of the bump. This is computed by taking
        // the derivative of the Bernoulli equation and finding the height
        // at the location of the maximum bump height where dz/dx = 0.
        //
        let h_middle = (q_in * q_in / g).powf(1.0 / 3.0);
        //
        // Iteration from the location of the maximum bump height to the beginning.
        //
        let mut h_water = vec![0.0; self.n_cells + 1];
        h_water[self.n_max] = h_middle;

        for n in (0..self.n_max).rev() {
            //
            // The coefficients of the Bernoulli cubic equation to solve
            // h^3 + a*h^2 + b*h + c = 0.
            //
            let coeffs =
                shallow_water::bernoulli_coefficients(q_in, h_middle, self.z[n], self.z_max, g);
            //
            // The bump height solutions.
            //
            let height = shallow_water::solve_cubic_equation(&coeffs);
            //
            // Find the solution which makes sense for this case.
            //
            h_water[n] = shallow_water::find_solution(&height, h_water[n + 1] * (1.0 + epsilon));
        }
        //
        // As the bump is symmetric and the values do not coincide with the maximum
        // of the topology, we make the solution decrease.
        //
        let coeffs = shallow_water::bernoulli_coefficients(
            q_in,
            h_middle,
            self.z[self.n_max + 1],
            self.z_max,
            g,
        );

        let height = shallow_water::solve_cubic_equation(&coeffs);
        h_water[self.n_max + 1] = shallow_water::find_solution(&height, h_middle);

        for n in (self.n_max + 2)..=self.n_cells {
            //
            // The coefficients of the Bernoulli cubic equation to solve
            // h^3 + a*h^2 + b*h + c = 0.
            //
            let coeffs =
                shallow_water::bernoulli_coefficients(q_in, h_middle, self.z[n], self.z_max, g);
            //
            // The bump height solutions.
            //
            let height = shallow_water::solve_cubic_equation(&coeffs);
            //
            // Find the solution which makes sense for this case.
            //
            h_water[n] = shallow_water::find_solution(&height, h_water[n - 1] * (1.0 - epsilon));
        }
        //
        // Output the results to the screen and a file.
        //
        if screen_output {
            println!("Maximum bump height                 = {:.6}", self.z_max);
            println!(
                "Location of maximum bump height     = {:.6}",
                self.x[self.n_max]
            );
            println!(
                "Water height at maximum bump height = {:.6}",
                h_water[self.n_max]
            );
            println!("==============================================");
            println!("x (m)     Z (m)     height (m)");
        }
        //
        // Use zip to iterate over corresponding elements in both vectors
        // and update the water height.
        //
        h_water.iter_mut().zip(&self.z).for_each(|(h, &z)| {
            *h += z;
        });
        //
        // Output results
        //
        for n in 1..=self.n_cells {
            //
            // Output the results to the screen.
            //
            if screen_output {
                println!("{:.6}  {:.6}  {:.6}", self.x[n], self.z[n], h_water[n]);
            }
            //
            // Output the results to a file.
            //
            writeln!(
                transcritical_no_shock_file,
                "{:.6}   {:.6}",
                self.x[n], h_water[n]
            );
        }

        transcritical_no_shock_file.flush().ok();
    }

    ///
    /// Computes the transcritical shock case and outputs the results.
    ///
    /// - Arguments:
    ///   - `q_in`: the inflow discharge [m²/s].
    ///   - `h_out`: the outflow height [m].
    ///   - `g`: the acceleration due to gravity [m/s²].
    ///   - `filename`: the name of the file to output the results.
    ///   - `screen_output`: whether to output the results to the screen.
    ///
    pub fn compute_transcritical_shock_case(
        &self,
        q_in: f64,
        h_out: f64,
        g: f64,
        filename: &str,
        screen_output: bool,
    ) {
        //
        // Trancritical shock flow over a bump.
        //
        if screen_output {
            println!();
            println!("Transcritical Shock Flow Over a Bump");
            println!("==============================================");
            println!("x (m)     Z (m)     height (m)");
        }
        //
        // Open a file to save the height of water.
        //
        let test_file = std::fs::File::create(filename);
        let mut transcritical_shock_file = match test_file {
            Ok(file) => file,
            Err(_) => {
                println!("Bump - Error opening file {}", filename);
                return;
            }
        };
        writeln!(transcritical_shock_file, "# x (m)    height (m)");
        //
        // Discretisation
        //
        let epsilon = 1.0 / self.n_cells as f64;
        let epsi = 10.0 / self.n_cells as f64;
        //
        // Water height at the top of the bump. This is computed by taking
        // the derivative of the Bernoulli equation and finding the height
        // at the location of the maximum bump height where dz/dx = 0.
        //
        let h_middle = (q_in * q_in / g).powf(1.0 / 3.0);
        //
        // Search for the shock location.
        //
        let mut rh_test = 100.0;
        let mut h_plus = 0.0;
        let mut h_minus = 0.0;
        let mut n_shock = self.n_max + 1;

        while rh_test > epsi && n_shock < self.n_cells {
            //
            // The coefficients of the Bernoulli cubic equation to solve
            // a*h^3 + b*h^2 + c*h + d = 0.
            //
            let coeffs_1 = shallow_water::bernoulli_coefficients(
                q_in,
                h_out,
                self.z[n_shock],
                self.z[self.n_cells],
                g,
            );
            //
            // Compute the h_plus water height.
            //
            let hplus = shallow_water::solve_cubic_equation(&coeffs_1);
            h_plus = shallow_water::find_solution(&hplus, h_out);
            if h_plus <= 1.0e-5 {
                n_shock += 1;
                continue;
            }
            //
            // The coefficients of the Bernoulli cubic equation to solve
            // a*h^3 + b*h^2 + c*h + d = 0.
            //
            let coeffs_2 = shallow_water::bernoulli_coefficients(
                q_in,
                h_middle,
                self.z[n_shock],
                self.z_max,
                g,
            );
            //
            // Compute the h_minus water height.
            //
            let hminus = shallow_water::solve_cubic_equation(&coeffs_2);
            h_minus = shallow_water::find_solution(&hminus, h_middle);
            if h_minus <= 1.0e-5 {
                n_shock += 1;
                continue;
            }
            //
            // Compute the Rankine-Hugoniot condition.
            //
            rh_test = Bump::rankine_hugoniot(h_plus, h_minus, q_in, g);
            n_shock += 1;
        }
        //
        // Water height with the shock location found.
        //
        let mut h_water = vec![0.0; self.n_cells + 1];
        h_water[n_shock] = h_minus;
        h_water[self.n_cells] = h_out;
        //
        // Iteration to compute the water height from the shock location to the
        // location of the inlet.
        //
        for n in (0..n_shock).rev() {
            //
            // The coefficients of the Bernoulli cubic equation to solve
            // a*h^3 + b*h^2 + c*h + d = 0.
            //
            let coeffs =
                shallow_water::bernoulli_coefficients(q_in, h_middle, self.z[n], self.z_max, g);
            //
            // The bump height solutions.
            //
            let height = shallow_water::solve_cubic_equation(&coeffs);
            //
            // Find the solution which makes sense for this case.
            //
            h_water[n] = shallow_water::find_solution(&height, h_water[n + 1] + epsilon);
        }
        //
        // Iteration to compute the water height from the outlet to the
        // location of the shock.
        //
        for n in (n_shock + 1..self.n_cells).rev() {
            //
            // The coefficients of the Bernoulli cubic equation to solve
            // a*h^3 + b*h^2 + c*h + d = 0.
            //
            let coeffs = shallow_water::bernoulli_coefficients(
                q_in,
                h_out,
                self.z[n],
                self.z[self.n_cells],
                g,
            );
            //
            // The bump height solutions.
            //
            let height = shallow_water::solve_cubic_equation(&coeffs);
            //
            // Find the solution which makes sense for this case.
            //
            h_water[n] = shallow_water::find_solution(&height, h_water[n + 1]);
        }
        //
        // Output the results to the screen and a file.
        //
        if screen_output {
            println!("Maximum bump height                 = {:.6}", self.z_max);
            println!(
                "Location of maximum bump height     = {:.6}",
                self.x[self.n_max]
            );
            println!(
                "Water height at maximum bump height = {:.6}",
                h_water[self.n_max]
            );
            println!(
                "Shock location (n, x)               = ({}, {})",
                n_shock, self.x[n_shock]
            );
            println!("Shock height minus                  = {:.6}", h_minus);
            println!("Shock height plus                   = {:.6}", h_plus);
            println!("==============================================");
            println!("x (m)     Z (m)     height (m)");
        }
        //
        // Use zip to iterate over corresponding elements in both vectors
        // and update the water height.
        //
        h_water.iter_mut().zip(&self.z).for_each(|(h, &z)| {
            *h += z;
        });
        //
        // Output results
        //
        for n in 1..=self.n_cells {
            //
            // Output the results to the screen.
            //
            if screen_output {
                println!("{:.6}  {:.6}  {:.6}", self.x[n], self.z[n], h_water[n]);
            }
            //
            // Output the results to a file.
            //
            writeln!(
                transcritical_shock_file,
                "{:.6}   {:.6}",
                self.x[n], h_water[n]
            );
        }

        transcritical_shock_file.flush().ok();
    }
}
