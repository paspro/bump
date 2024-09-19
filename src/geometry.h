/**
 * @file geometry.h
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

#include <algorithm>

/**
 * @brief The geometry of the bump.
 *
 * @param x The position.
 *
 * @return The height of the bump at position x.
 */
inline double
bump_geometry(double x)
{
    return std::max(0.0, 0.2 - 0.05 * (x - 10.0) * (x - 10.0));
}

/**
 * @brief The maximum height of the bump.
 *
 * @return The maximum height of the bump.
 */
inline double
maximum_bump_height()
{
    return 0.2;
}
