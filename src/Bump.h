/*!
 * @file Bump.h
 * @author Panos Asproulis <p.asproulis@icloud.com>
 * @version 2.1
 *
 * @brief Computes the analytical solution of the shallow water flow over a
 * bump.
 */

#pragma once

#include <algorithm>
#include <string>
#include <valarray>

/*!
 * @class Bump
 * @brief A class to compute the analytical solution of the shallow water flow
 * over a bump.
 *
 * The Bump class provides static methods to compute the geometry of a bump
 * and the height of water for the subcritical and supercritical cases. This
 * class cannot be instantiated or copied.
 */
class Bump
{

  public:
    /*!
     * @brief Deleted default constructor.
     */
    Bump() = delete;
    /*!
     * @brief Deleted destructor.
     */
    ~Bump() = delete;
    /*!
     * @brief Deleted copy constructor.
     */
    Bump(const Bump&) = delete;
    /*!
     * @brief Deleted copy assignment operator.
     */
    Bump& operator=(const Bump&) = delete;

    /*!
     * @brief Initialise the geometry of the bump.
     *
     * @param length    Length of the domain [m].
     * @param n_cells   Number of discretisation cells.
     */
    static void initialise_geometry(double length, std::size_t n_cells);

    /*!
     * @brief Computes the Rankine-Hugoniot condition for shock conditions.
     *
     * @param h_plus    Height of water on right side of shock [m].
     * @param h_minus   Height of water on left side of shock [m].
     * @param q         Flow rate of water [m²/s].
     * @param g         Gravitational acceleration [m/s²].
     *
     * @return Value of Rankine-Hugoniot condition (should be zero if satisfied).
     */
    inline static double rankine_hugoniot(double h_plus, double h_minus, double q, double g)
    {
        return std::fabs(q * q * (1.0 / h_plus - 1.0 / h_minus) +
                         0.5 * g * (h_plus - h_minus) * (h_plus + h_minus));
    }

    /*!
     * @brief Compute the subcritical case and outputs the results.
     *
     * @param q_in           Water flow rate at inlet [m²/s].
     * @param h_out          Water height at outlet [m].
     * @param g              Gravitational acceleration [m/s²].
     * @param filename       Name of the file to save results.
     * @param screen_output  Flag to output results to screen.
     */
    static void compute_subcritical_case(double q_in,
                                         double h_out,
                                         double g,
                                         const std::string& filename,
                                         bool screen_output = false);

    /*!
     * @brief Compute the transcritical no-shock case and outputs the results.
     *
     * @param q_in           Water flow rate at inlet [m²/s].
     * @param g              Gravitational acceleration [m/s²].
     * @param filename       Name of the file to save results.
     * @param screen_output  Flag to output results to screen.
     */
    static void compute_transcritical_no_shock_case(double q_in,
                                                    double g,
                                                    const std::string& filename,
                                                    bool screen_output = false);

    /*!
     * @brief Compute the transcritical shock case and outputs the results.
     *
     * @param q_in           Water flow rate at inlet [m²/s].
     * @param h_out          Water height at outlet [m].
     * @param g              Gravitational acceleration [m/s²].
     * @param filename       Name of the file to save results.
     * @param screen_output  Flag to output results to screen.
     */
    static void compute_transcritical_shock_case(double q_in,
                                                 double h_out,
                                                 double g,
                                                 const std::string& filename,
                                                 bool screen_output = false);

  private:
    /*!
     * @brief Computes the geometry of the bump.
     *
     * @param x  Position along domain [m].
     *
     * @return Height of bump at position x [m].
     */
    inline static double bump_geometry(double x)
    {
        return std::max(0.0, 0.2 - 0.05 * (x - 10.0) * (x - 10.0));
    }

  private:
    /*!
     * @brief The length of the domain.
     */
    static double length; ///< Domain length [m]
    /*!
     * @brief The number of cells to use for the discretisation process.
     */
    static std::size_t n_cells; ///< Number of discretisation cells
    /*!
     * @brief The x-coordinates of the cell centers.
     */
    static std::valarray<double> x; ///< Cell center x-coordinates [m]
    /*!
     * @brief The height of the bump at the cell centers.
     */
    static std::valarray<double> z; ///< Bump height at cell centers [m]
    /*!
     * @brief The maximum height of the bump.
     */
    static double z_max; ///< Maximum bump height [m]
    /*!
     * @brief The location of the maximum height of the bump.
     */
    static std::size_t n_max; ///< Index of maximum bump height location

}; // class Bump
