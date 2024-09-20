/**
 * @file driver.cpp
 * @author Panos Asproulis <p.asproulis@icloud.com>
 * @version 1.0
 * @date 2024-09-19
 *
 * @brief Computes the analytical solution of the shallow water flow over a
 * bump. This program calculates the height of water for the case of shallow
 * flow over a bump. It handles both subcritical and transcritical flow cases.
 *
 * @copyright (c) Renuda (UK) Ltd., 2024
 */

#include <iostream>

#include "Bump.h"

int
main(void)
{
    std::size_t n_cells;
    std::string filename;
    double q_in, h_out, length;
    //
    // Title
    //
    std::cout << std::endl
              << "Bump - Computation of the height of water for the case of\n"
                 "shallow flow over a bump."
              << std::endl
              << std::endl;

    length  = 25.0; // length of the domain (m)
    n_cells = 200;  // number of cells for the computation
    Bump::initialise_geometry(length, n_cells);
    
    //
    // Subcritical flow over a bump
    //
    q_in     = 4.42;                   // water flow rate at x = 0 (m^2/s)
    h_out    = 2.00;                   // water height at the outlet (m)
    filename = "subcritical_bump.dat"; // name of the file to save the results
    Bump::compute_subcritical_case(
        q_in, h_out, length, 9.81, n_cells, filename, true);

    //
    // Transcritical no-shock flow over a bump
    //
    q_in     = 1.53; // water flow rate at x = 0 (m^2/s)
    filename = "transcritical_no_shock_bump.dat"; // name of the file to save
                                                  // the results
    Bump::compute_transcritical_no_shock_case(
        q_in, length, 9.81, n_cells, filename, true);

    //
    // Transcritical shock flow over a bump
    //
    q_in     = 0.18; // water flow rate at x = 0 (m^2/s)
    h_out    = 0.33; // water height at the outlet (m)
    filename = "transcritical_shock_bump.dat"; // name of the file to save
                                               // the results
    Bump::compute_transcritical_shock_case(
        q_in, h_out, length, 9.81, n_cells, filename, true);

    return (0);
}
