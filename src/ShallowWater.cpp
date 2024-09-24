/**
 * @file ShallowWater.cpp
 * @author Panos Asproulis <p.asproulis@icloud.com>
 * @version 2.1
 * @date 2024-09-24
 *
 * @brief Computes the analytical solution of the shallow water flow over a
 * bump.
 */

#include "ShallowWater.h"

//
// Compute the Bernoulli coefficients i.e. the coefficient of the cubic form
// of the Bernoulli equation a*h^3 + b*h^2 + c*h + d = 0 for the shallow
// water flow between an inflow and an outflow location.
//
// @param q_in The inflow water rate.
// @param h_out The water height at the outflow.
// @param z_in The elevation at the inflow.
// @param z_out The elevation at the outflow.
// @param g Gravity.
//
// @return The coefficients of the cubic form of the Bernoulli equation.
//
std::tuple<double, double, double, double>
ShallowWater::bernoulli_coefficients(
    double q_in, double h_out, double z_in, double z_out, double g)
{
    const double a = 1.0;
    const double b
        = -(q_in * q_in / (2.0 * g * h_out * h_out) + h_out - (z_in - z_out));
    const double c = 0.0;
    const double d = q_in * q_in / (2.0 * g);

    return std::make_tuple(a, b, c, d);
}

/**
 * @brief Solve a cubic equation of the form: a*h^3 + b*h^2 + c*h + d = 0
 * using Cardano's method. The cubic equation has three roots.
 *
 * @param coefficients The coefficient of the cubic equation.
 *
 * @return A vector of complex numbers representing the roots of the
 * cubic equation.
 */
std::vector<std::complex<double>>
ShallowWater::solve_cubic_equation(std::tuple<double, double, double, double> coefficients)
{
    double a = std::get<0>(coefficients);
    double b = std::get<1>(coefficients);
    double c = std::get<2>(coefficients);
    double d = std::get<3>(coefficients);
    //
    // Normalize the coefficients
    //
    b /= a;
    c /= a;
    d /= a;
    //
    // Calculate the discriminant
    //
    const double delta0 = b * b - 3.0 * c;
    const double delta1 = 2.0 * b * b * b - 9.0 * b * c + 27.0 * d;

    std::complex<double> C
        = std::pow((delta1
                    + std::sqrt(std::complex<double>(
                        delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)))
                       / 2.0,
                   1.0 / 3.0);

    const std::complex<double> omega(-0.5,
                                     std::sqrt(3.0)
                                         / 2.0); // Cube roots of unity
    //
    // Compute the roots
    //
    std::vector<std::complex<double>> roots(3);
    for (size_t k = 0; k < 3; ++k)
    {
        roots[k] = -1.0 / 3.0
                   * (b + std::pow(omega, k) * C
                      + delta0 / (std::pow(omega, k) * C));
    }

    return roots;
}

/**
 * @brief Finds the solution which makes sense for a particular case given the
 * three roots of the cubic equation.
 *
 * @param roots The three roots of the cubic equation
 * @param h_near The height of water at the previously computed location
 *
 * @return The height of water
 */
double
ShallowWater::find_solution(const std::vector<std::complex<double>>& roots, double h_near)
{
    //
    // Ignore the roots with imaginary parts and negative
    // real parts
    //
    std::vector<double> real_roots;
    for (const auto& root : roots)
    {
        if (std::fabs(root.imag()) <= 1.0e-6 && root.real() >= 0.0)
        {
            real_roots.push_back(root.real());
        }
    }
    //
    // Select the root which is closest to the previously
    // computed height of water
    //
    double solution = 0.0;
    double min_diff = std::numeric_limits<double>::max();

    for (const auto& root : real_roots)
    {
        const double diff = std::fabs(root - h_near);
        if (diff < min_diff)
        {
            solution = root;
            min_diff = diff;
        }
    }

    return solution;
}
