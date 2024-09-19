/**
 * @file bernoulli.h
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

#include <tuple>

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
inline std::tuple<double, double, double, double>
bernoulli_coefficients(
    double q_in, double h_out, double z_in, double z_out, double g)
{
    const double a = 1.0;
    const double b
        = -(q_in * q_in / (2.0 * g * h_out * h_out) + h_out - (z_in - z_out));
    const double c = 0.0;
    const double d = q_in * q_in / (2.0 * g);

    return std::make_tuple(a, b, c, d);
}