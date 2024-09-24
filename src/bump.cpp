/**
 * @file Bump.cpp
 * @author Panos Asproulis <p.asproulis@icloud.com>
 * @version 2.1
 * @date 2024-09-24
 *
 * @brief Computes the analytical solution of the shallow water flow over a
 * bump.
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Bump.h"
#include "ShallowWater.h"

//
// Definition of the static variables
//
double Bump::length;
std::size_t Bump::n_cells;
std::valarray<double> Bump::x;
std::valarray<double> Bump::z;
double Bump::z_max;
std::size_t Bump::n_max;

/**
 * @brief Initialise the geometry of the bump.
 *
 * @param length Length of the domain.
 * @param n_cells Number of discretisation cells.
 */
void
Bump::initialise_geometry(double length, std::size_t n_cells)
{
    //
    // Discretisation of the domain
    //
    Bump::length  = length;
    Bump::n_cells = n_cells;

    const double dx = length / double(n_cells);
    Bump::x.resize(n_cells + 1);
    Bump::z.resize(n_cells + 1);

    for (int n = 0; n <= n_cells; n++)
    {
        Bump::x[n] = (double(n) - 0.5) * dx;
        Bump::z[n] = Bump::bump_geometry(Bump::x[n]);
    }
    //
    // Compute the maximum bump height and its location
    //
    Bump::z_max = 0.0;
    Bump::n_max = 0;

    for (int n = 0; n <= n_cells; n++)
    {
        if (Bump::z[n] > Bump::z_max)
        {
            Bump::z_max = Bump::z[n];
            Bump::n_max = n;
        }
    }
    Bump::z_max = 0.2;
}

/**
 * @brief Compute the subcritical case and outputs the results.
 *
 * @param q_in Water flow rate at x = 0.
 * @param h_out Water height at the outlet.
 * @param g Gravity.
 * @param filename Name of the file to save the results.
 * @param screen_output Flag to output the results to the screen.
 */
void
Bump::compute_subcritical_case(double q_in,
                               double h_out,
                               double g,
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
        std::cout << "x (m)     Z (m)     height (m)" << std::endl;
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
    // Iteration from the end of the field to the beginning
    // with application of the Bernoulli equation at each cell
    //
    std::valarray<double> h_water(0.0, n_cells + 1);
    h_water[n_cells] = h_out;

    for (int n = Bump::n_cells - 1; n >= 0; n--)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs = ShallowWater::bernoulli_coefficients(
            q_in, h_out, Bump::z[n], 0.0, g);
        //
        // The bump height solutions
        //
        const auto height = ShallowWater::solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = ShallowWater::find_solution(height, h_water[n + 1]);
    }
    for (int n = 0; n <= Bump::n_cells; n++)
    {
        h_water[n] += Bump::z[n];
    }
    //
    // Output the results to the screen and a file
    //
    for (int n = 1; n <= Bump::n_cells; n++)
    {
        //
        // Output the results to the screen
        //
        if (screen_output)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << Bump::x[n] << "  " << z[n] << "  " << h_water[n] << std::endl;
        }
        //
        // Output the results to a file
        //
        if (!filename.empty())
        {
            subcritical_file << std::fixed << std::setprecision(6);
            subcritical_file << Bump::x[n] << "   " << h_water[n] << std::endl;
        }
    }
    //
    // Close the file
    //
    if (!filename.empty())
    {
        subcritical_file.flush();
        subcritical_file.close();
    }
}

/**
 * @brief Compute the transcritical no-shock case and outputs the results.
 *
 * @param q_in Water flow rate at x = 0.
 * @param g Gravity.
 * @param filename Name of the file to save the results.
 * @param screen_output Flag to output the results to the screen.
 */
void
Bump::compute_transcritical_no_shock_case(double q_in,
                                          double g,
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
    const double epsilon = 1.0 / double(Bump::n_cells);
    //
    // Water height at the top of the bump. This is computed by taking
    // the derivative of the Bernoulli equation and finding the height
    // at the location of the maximum bump height where dz/dx = 0.
    //
    const double h_middle = pow(q_in * q_in / g, 1.0 / 3.0);
    //
    // Iteration from the location of the maximum bump height to the beginning
    //
    std::valarray<double> h_water(0.0, Bump::n_cells + 1);
    h_water[Bump::n_max] = h_middle;

    for (int n = Bump::n_max - 1; n >= 0; n--)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // h^3 + a*h^2 + b*h + c = 0
        //
        const auto coeffs = ShallowWater::bernoulli_coefficients(
            q_in, h_middle, Bump::z[n], Bump::z_max, g);
        //
        // The bump height solutions
        //
        auto height = ShallowWater::solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = ShallowWater::find_solution(
            height, h_water[n + 1] * (1.0 + epsilon));
    }
    //
    // As the bump is symmetric and the values do not coincide with the maximum
    // of the topology, we make the solution deacrease.
    //
    const auto coeffs = ShallowWater::bernoulli_coefficients(
        q_in, h_middle, Bump::z[Bump::n_max + 1], Bump::z_max, g);

    auto height              = ShallowWater::solve_cubic_equation(coeffs);
    h_water[Bump::n_max + 1] = ShallowWater::find_solution(height, h_middle);

    for (int n = Bump::n_max + 2; n <= Bump::n_cells; n++)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // h^3 + a*h^2 + b*h + c = 0
        //
        const auto coeffs = ShallowWater::bernoulli_coefficients(
            q_in, h_middle, Bump::z[n], Bump::z_max, g);
        //
        // The bump height solutions
        //
        auto height = ShallowWater::solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = ShallowWater::find_solution(
            height, h_water[n + 1] * (1.0 - epsilon));
    }
    //
    // Output the results to the screen and a file
    //
    if (screen_output)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Maximum bump height                 = " << Bump::z_max
                  << std::endl;
        std::cout << "Location of maximum bump height     = " << Bump::x[Bump::n_max]
                  << std::endl;
        std::cout << "Water height at maximum bump height = " << h_water[Bump::n_max]
                  << std::endl;
        std::cout << "=============================================="
                  << std::endl;
        std::cout << "x (m)     Z (m)     height (m)" << std::endl;
    }
    for (int n = 0; n <= Bump::n_cells; n++)
    {
        h_water[n] += Bump::z[n];
    }
    for (int n = 1; n <= n_cells; n++)
    {
        //
        // Output the results to the screen
        //
        if (screen_output)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << Bump::x[n] << "  " << z[n] << "  " << h_water[n] << std::endl;
        }
        //
        // Output the results to a file
        //
        if (!filename.empty())
        {
            transcritical_no_shock_file << std::fixed << std::setprecision(6);
            transcritical_no_shock_file << Bump::x[n] << "   " << h_water[n]
                                        << std::endl;
        }
    }
    //
    // Close the file
    //
    if (!filename.empty())
    {
        transcritical_no_shock_file.flush();
        transcritical_no_shock_file.close();
    }
}

/**
 * @brief Compute the transcritical shock case and outputs the results.
 *
 * @param q_in Water flow rate at x = 0.
 * @param h_out Water height at the outlet.
 * @param g Gravity.
 * @param filename Name of the file to save the results.
 * @param screen_output Flag to output the results to the screen.
 */
void
Bump::compute_transcritical_shock_case(double q_in,
                                       double h_out,
                                       double g,
                                       std::string filename,
                                       bool screen_output)
{
    //
    // Transcritical No-Shock flow over a bump
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
    const double epsilon = 1.0 / double(Bump::n_cells);
    const double epsi    = 10.0 / double(Bump::n_cells);
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
    std::size_t n_shock = n_max + 1;

    while (rh_test > epsi && n_shock < Bump::n_cells)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs1 = ShallowWater::bernoulli_coefficients(
            q_in, h_out, Bump::z[n_shock], Bump::z[n_cells], g);
        //
        // Compute the h_plus water height
        //
        const auto hplus = ShallowWater::solve_cubic_equation(coeffs1);
        h_plus           = ShallowWater::find_solution(hplus, h_out);
        if (h_plus <= 1.e-5) {
            n_shock++;
            continue;
        }
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs2 = ShallowWater::bernoulli_coefficients(
            q_in, h_middle, Bump::z[n_shock], Bump::z_max, g);
        //
        // Compute the h_minus water height
        //
        const auto hminus = ShallowWater::solve_cubic_equation(coeffs2);
        h_minus           = ShallowWater::find_solution(hminus, h_middle);
        if (h_minus <= 1.e-5)
        {
            n_shock++;
            continue;
        };
        //
        // Compute the Rankine-Hugoniot condition
        //
        rh_test = Bump::rankine_hugoniot(h_plus, h_minus, q_in, g);
        n_shock++;
    }
    //
    // Water height with the shock location found
    //
    std::valarray<double> h_water(0.0, Bump::n_cells + 1);
    h_water[n_shock] = h_minus;
    h_water[Bump::n_cells] = h_out;
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
        const auto coeffs = ShallowWater::bernoulli_coefficients(
            q_in, h_middle, Bump::z[n], Bump::z_max, g);
        //
        // The bump height solutions
        //
        const auto height = ShallowWater::solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n]
            = ShallowWater::find_solution(height, h_water[n + 1] + epsilon);
    }
    //
    // Iteration to compute the water height from the outlet to the
    // location of the shock
    //
    for (int n = Bump::n_cells - 1; n > n_shock; n--)
    {
        //
        // The coefficients of the Bernoulli cubic equation to solve
        // a*h^3 + b*h^2 + c*h + d = 0
        //
        const auto coeffs = ShallowWater::bernoulli_coefficients(
            q_in, h_out, Bump::z[n], Bump::z[Bump::n_cells], g);
        //
        // The bump height solutions
        //
        const auto height = ShallowWater::solve_cubic_equation(coeffs);
        //
        // Find the solution which makes sense for this case
        //
        h_water[n] = ShallowWater::find_solution(height, h_water[n + 1]);
    }
    //
    // Output the results to the screen and a file
    //
    if (screen_output)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Maximum bump height                 = " << Bump::z_max
                  << std::endl;
        std::cout << "Location of maximum bump height     = " << Bump::x[Bump::n_max]
                  << std::endl;
        std::cout << "Water height at maximum bump height = " << h_water[Bump::n_max]
                  << std::endl;
        std::cout << "Shock location (n, x)               = "
                  << "(" << n_shock << ", " << Bump::x[n_shock] << ")" << std::endl;
        std::cout << "Shock height minus                  = " << h_minus
                  << std::endl;
        std::cout << "Shock height plus                   = " << h_plus
                  << std::endl;
        std::cout << "=============================================="
                  << std::endl;
        std::cout << "x (m)     Z (m)     height (m)" << std::endl;
    }
    for (int n = 0; n <= Bump::n_cells; n++)
    {
        h_water[n] += Bump::z[n];
    }
    for (int n = 1; n <= Bump::n_cells; n++)
    {
        //
        // Output the results to the screen
        //
        if (screen_output)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << Bump::x[n] << "  " << z[n] << "  " << h_water[n] << std::endl;
        }
        //
        // Output the results to a file
        //
        if (!filename.empty())
        {
            transcritical_shock_file << std::fixed << std::setprecision(6);
            transcritical_shock_file << Bump::x[n] << "   " << h_water[n]
                                     << std::endl;
        }
    }
    //
    // Close the file
    //
    if (!filename.empty())
    {
        transcritical_shock_file.flush();
        transcritical_shock_file.close();
    }
}
