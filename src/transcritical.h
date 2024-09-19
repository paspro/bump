/**
 * @file transcritical.h
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

#include <cstdlib>
#include <string>

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
                                    bool screen_output = false);

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
                                 bool screen_output = false);