/**
 * @file ShallowWater.h
 * @author Panos Asproulis <p.asproulis@icloud.com>
 * @version 2.0
 * @date 2024-09-20
 *
 * @brief Computes the analytical solution of the shallow water flow over a
 * bump.
 */

#pragma once

#include <complex>
#include <tuple>
#include <vector>

class ShallowWater
{

public:
    /**
     * @brief Deleted default constructor
     */
    ShallowWater() = delete;
    /**
     * @brief Deleted destructor
     */
    ~ShallowWater() = delete;
    /**
     * @brief Deleted copy constructor
     */
    ShallowWater(const ShallowWater&) = delete;
    /**
     * @brief Deleted copy assignment operator.
     */
    ShallowWater&
    operator=(const ShallowWater&) = delete;

public:
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
    static std::tuple<double, double, double, double>
    bernoulli_coefficients(
        double q_in, double h_out, double z_in, double z_out, double g);

    /**
     * @brief Solve a cubic equation of the form: a*h^3 + b*h^2 + c*h + d = 0
     * using Cardano's method. The cubic equation has three roots.
     *
     * @param coefficients The coefficient of the cubic equation.
     *
     * @return A vector of complex numbers representing the roots of the
     * cubic equation.
     */
    static std::vector<std::complex<double>>
    solve_cubic_equation(
        std::tuple<double, double, double, double> coefficients);

    /**
     * @brief Find the solution which makes sense for a particular case
     * given the three roots of the cubic equation. The function returns
     * the height of water.
     *
     * @param roots The three roots of the cubic equation.
     * @param h_near The height of water at the previously computed location.
     *
     * @return The height of water.
     */
    static double
    find_solution(const std::vector<std::complex<double>>& roots,
                  double h_near);

}; // class ShallowWater
