/**
 * @file transcritical.cpp
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
#include "transcritical.h"

/**
 * @brief It computes the Rankine-Hugoniot condition and returns the value of
 * the equation. If the input is correct then the value should be zero.
 *
 * @param h_plus The height of the water on the right side of the shock.
 * @param h_minus The height of the water on the left side of the shock.
 * @param q The flow rate of the water.
 * @param g Gravity constant
 *
 * @return The value of the Rankine-Hugoniot condition.
 */
inline double
rankine_hugoniot(double h_plus, double h_minus, double q, double g)
{
    return abs(q * q * (1.0 / h_plus - 1.0 / h_minus)
               + 0.5 * g * (h_plus - h_minus) * (h_plus + h_minus));
}

/**
 * @brief Compute the transcritical no-shock case and outputs the results.
 *
 * @param q_in Water flow rate at x = 0.
 * @param length Length of the domain.
 * @param g Gravity.
 * @param n_cells Number of discretisation cells.
 * @param filename Name of the file to save the results.
 * @param screen_output Flag to output the results to the screen.
 */
void
compute_transcritical_no_shock_case(double q_in,
                                    double length,
                                    double g,
                                    size_t n_cells,
                                    std::string filename,
                                    bool screen_output)
{
    //
    // Transcritical No-Shock flow over a bump
    //
    if (screen_output)
    {
        std::cout << std::endl;
        std::cout << "Transcritical No-Shock Flow Over a Bump" << std::endl;
        std::cout << "=============================================="
                  << std::endl;
    }
    //
    // Open a file to save the height of water
    //
    std::ofstream transcritical_no_shock_file;
    if (!filename.empty())
    {
        transcritical_no_shock_file.open(filename, std::ios::out);
        transcritical_no_shock_file << "# x (m)    height (m)" << std::endl;
    }
    //
    // Discretisation
    //
    const double epsilon = 1.0 / double(n_cells);
    const double dx      = length / double(n_cells);
    std::valarray<double> x(0.0, n_cells + 1);
    std::valarray<double> z(0.0, n_cells + 1);

    for (int n = 0; n <= n_cells; n++)
    {
        x[n] = (double(n) - 0.5) * dx;
        z[n] = bump_geometry(x[n]);
    }
    //
    // Find the value and the location of the maximum bump height
    //
    double z_max      = 0.0;
    std::size_t n_max = 0;

    for (int n = 0; n <= n_cells; n++)
    {
        if (z[n] > z_max)
        {
            z_max = z[n];
            n_max = n;
        }
    }
    z_max = maximum_bump_height();
    //
    // Water height at the top of the bump. This is computed by taking
    // the derivative of the Bernoulli equation and finding the height
    // at the location of the maximum bump height where dz/dx = 0.
    //
    const double h_middle = pow(q_in * q_in / g, 1.0 / 3.0);
    //
    // Iteration from the location of the maximum bump height to the beginning
    //
    std::valarray<double> h_water(0.0, n_cells + 1);
    h_water[n_max] = h_middle;

    for (int n = n_max - 1; n >= 0; n--)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // h^3 + a*h^2 + b*h + c = 0
        //
        const auto coeffs
            = bernoulli_coefficients(q_in, h_middle, z[n], z_max, g);
        //
        // The bump height solutions
        //
        auto height = solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = find_solution(height, h_water[n + 1] * (1.0 + epsilon));
    }
    //
    // As the bump is symmetric and the values do not coincide with the maximum
    // of the topology, we make the solution deacrease.
    //
    const auto coeffs
        = bernoulli_coefficients(q_in, h_middle, z[n_max + 1], z_max, g);
    auto height        = solve_cubic_equation(coeffs);
    h_water[n_max + 1] = find_solution(height, h_middle);

    for (int n = n_max + 2; n <= n_cells; n++)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // h^3 + a*h^2 + b*h + c = 0
        //
        const auto coeffs
            = bernoulli_coefficients(q_in, h_middle, z[n], z_max, g);
        //
        // The bump height solutions
        //
        auto height = solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = find_solution(height, h_water[n + 1] * (1.0 - epsilon));
    }
    //
    // Output the results to the screen and a file
    //
    if (screen_output)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Maximum bump height                 = " << z_max
                  << std::endl;
        std::cout << "Location of maximum bump height     = " << x[n_max]
                  << std::endl;
        std::cout << "Water height at maximum bump height = " << h_water[n_max]
                  << std::endl;
        std::cout << "=============================================="
                  << std::endl;
        std::cout << "x (m)      height (m)" << std::endl;
    }
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
            transcritical_no_shock_file << std::fixed << std::setprecision(6);
            transcritical_no_shock_file << x[n] << "   " << h_water[n]
                                        << std::endl;
        }
    }
    //
    // Close the file
    //
    if (!filename.empty())
    {
        transcritical_no_shock_file.close();
    }
}

/**
 * @brief Compute the transcritical shock case and outputs the results.
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
compute_transcritical_shock_case(double q_in,
                                 double h_out,
                                 double length,
                                 double g,
                                 size_t n_cells,
                                 std::string filename,
                                 bool screen_output)
{
    //
    // Transcritical Shock flow over a bump
    //
    if (screen_output)
    {
        std::cout << std::endl;
        std::cout << "Transcritical Shock Flow Over a Bump" << std::endl;
        std::cout << "=============================================="
                  << std::endl;
    }
    //
    // Open a file to save the height of water
    //
    std::ofstream transcritical_shock_file;
    if (!filename.empty())
    {
        transcritical_shock_file.open(filename, std::ios::out);
        transcritical_shock_file << "# x (m)    height (m)" << std::endl;
    }
    //
    // Discretisation
    //
    const double dx      = length / double(n_cells);
    const double epsilon = 1.0 / double(n_cells);
    const double epsi    = 10.0 / double(n_cells);
    std::valarray<double> x(0.0, n_cells + 1);
    std::valarray<double> z(0.0, n_cells + 1);

    for (int n = 0; n <= n_cells; n++)
    {
        x[n] = (double(n) - 0.5) * dx;
        z[n] = bump_geometry(x[n]);
    }
    //
    // Find the location and value of the maximum bump height
    //
    double z_max      = 0.0;
    std::size_t n_max = 0;

    for (int n = 0; n <= n_cells; n++)
    {
        if (z[n] > z_max)
        {
            z_max = z[n];
            n_max = n;
        }
    }
    z_max = maximum_bump_height();
    //
    // Water height at the top of the bump. This is computed by taking
    // the derivative of the Bernoulli equation and finding the height
    // at the location of the maximum bump height where dz/dx = 0.
    //
    const double h_middle = pow(q_in * q_in / g, 1.0 / 3.0);
    //
    // Search for the shock location
    //
    double rh_test = 100.0;
    double h_plus = 0.0, h_minus = 0.0;
    std::size_t n_shock = n_max;

    while (rh_test > epsi && n_shock < n_cells)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs1
            = bernoulli_coefficients(q_in, h_out, z[n_shock], z[n_cells], g);
        //
        // Compute the h_plus water height
        //
        const auto hplus = solve_cubic_equation(coeffs1);
        h_plus           = find_solution(hplus, h_out);
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs2
            = bernoulli_coefficients(q_in, h_middle, z[n_shock], z_max, g);
        //
        // Compute the h_minus water height
        //
        const auto hminus = solve_cubic_equation(coeffs2);
        h_minus           = find_solution(hminus, h_middle);
        //
        // Compute the Rankine-Hugoniot condition
        //
        rh_test = rankine_hugoniot(h_plus, h_minus, q_in, g);
        n_shock++;
    }
    //
    // Water height with the shock location found
    //
    std::valarray<double> h_water(0.0, n_cells + 1);
    h_water[n_shock] = h_minus;
    h_water[n_cells] = h_out;
    //
    // Iteration to compute the water height from the shock location to the
    // location of the inlet.
    //
    for (int n = n_shock - 1; n >= 0; n--)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs
            = bernoulli_coefficients(q_in, h_middle, z[n], z_max, g);
        //
        // The bump height solutions
        //
        const auto height = solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = find_solution(height, h_water[n + 1] + epsilon);
    }
    //
    // Iteration to compute the water height from the outlet to the
    // location of the shock
    //
    for (int n = n_cells - 1; n > n_shock; n--)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs
            = bernoulli_coefficients(q_in, h_out, z[n], z[n_cells], g);
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
    if (screen_output)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Maximum bump height                 = " << z_max
                  << std::endl;
        std::cout << "Location of maximum bump height     = " << x[n_max]
                  << std::endl;
        std::cout << "Water height at maximum bump height = " << h_water[n_max]
                  << std::endl;
        std::cout << "Shock location                      = " << x[n_shock]
                  << std::endl;
        std::cout << "Shock height minus                  = " << h_minus
                  << std::endl;
        std::cout << "Shock height plus                   = " << h_plus
                  << std::endl;
        std::cout << "=============================================="
                  << std::endl;
        std::cout << "x (m)      height (m)" << std::endl;
    }
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
            transcritical_shock_file << std::fixed << std::setprecision(6);
            transcritical_shock_file << x[n] << "   " << h_water[n]
                                     << std::endl;
        }
    }
    //
    // Close the file
    //
    if (!filename.empty())
    {
        transcritical_shock_file.close();
    }
}