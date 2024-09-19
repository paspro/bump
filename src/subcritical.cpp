/**
 * @file subcritical.cpp
 * @author Panos Asproulis <p.asproulis@icloud.com>
 * @version 1.0
 * @date 2024-09-19
 *
 * @brief Computes the analytical solution of the shallow water flow over a
 * bump.
 *
 * @copyright (c) Renuda (UK) Ltd., 2024
 */

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <valarray>
#include <vector>

#include "bernoulli.h"
#include "cubic.h"
#include "geometry.h"
#include "subcritical.h"

/**
 * @brief Compute the subcritical case and outputs the results.
 *
 * @param q_in Water flow rate at x = 0.
 * @param h_out Water height at the outlet.
 * @param length Length of the domain.
 * @param g Gravity.
 * @param n_cells Number of discretisation cells.
 * @param filename Name of the file to save the results.
 * @param screen_output Flag to output the results to the screen.
 */
void
compute_subcritical_case(double q_in,
                         double h_out,
                         double length,
                         double g,
                         size_t n_cells,
                         std::string filename,
                         bool screen_output)
{
    //
    // Subcritical flow over a bump
    //
    if (screen_output)
    {
        std::cout << std::endl;
        std::cout << "Subcritical Flow Over a Bump" << std::endl;
        std::cout << "=============================================="
                  << std::endl;
        std::cout << "x (m)      height (m)" << std::endl;
    }
    //
    // Open a file to save the height of water
    //
    std::ofstream subcritical_file;
    if (!filename.empty())
    {
        subcritical_file.open(filename, std::ios::out);
        subcritical_file << "# x (m)    height (m)" << std::endl;
    }
    //
    // Discretisation
    //
    const double dx = length / double(n_cells);
    std::valarray<double> x(0.0, n_cells + 1);
    std::valarray<double> z(0.0, n_cells + 1);

    for (int n = 0; n <= n_cells; n++)
    {
        x[n] = (double(n) - 0.5) * dx;
        z[n] = bump_geometry(x[n]);
    }
    //
    // Iteration from the end of the field to the beginning
    //
    std::valarray<double> h_water(0.0, n_cells + 1);
    h_water[n_cells] = h_out;

    for (int n = n_cells - 1; n >= 0; n--)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs = bernoulli_coefficients(q_in, h_out, z[n], 0.0, g);
        //
        // The bump height solutions
        //
        const auto height = solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = find_solution(height, h_water[n + 1]);
    }
    //
    // Output the results to the screen and a file
    //
    for (int n = 1; n <= n_cells; n++)
    {
        //
        // Output the results to the screen
        //
        if (screen_output)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << x[n] << "   " << h_water[n] << std::endl;
        }
        //
        // Output the results to a file
        //
        if (!filename.empty())
        {
            subcritical_file << std::fixed << std::setprecision(6);
            subcritical_file << x[n] << "   " << h_water[n] << std::endl;
        }
    }
    //
    // Close the file
    //
    if (!filename.empty())
    {
        subcritical_file.close();
    }
}
