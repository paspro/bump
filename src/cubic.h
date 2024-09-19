/**
 * @file cubic.h
 * @author Panos Asproulis <p.asproulis@icloud.com>
 * @version 1.0
 * @date 2024-09-19
 *
 * @brief Computes the analytical solution of the shallow water flow over a
 * bump.
 *
 * @copyright (c) Renuda (UK) Ltd., 2024
 */

#pragma once

#include <complex>
#include <tuple>
#include <vector>

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
solve_cubic_equation(std::tuple<double, double, double, double> coefficients);

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
double
find_solution(const std::vector<std::complex<double>>& roots, double h_near);