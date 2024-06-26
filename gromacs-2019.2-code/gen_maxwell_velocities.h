/* Copyright (c) 2023, Huiyu Li and Ao Ma, University of Illinois at Chicago.
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef GMX_MAXWELL_VELOCITIES
#define GMX_MAXWELL_VELOCITIES

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;

/*! \brief
 * Generate Maxwellian velocities.
 *
 * \param[in] tempi Temperature to generate around
 * \param[in] seed  Random number generator seed
 * \param[in] mtop  Molecular Topology
 * \param[out] v    Velocities
 */
void maxwell_speed(real tempi, unsigned int seed,
                   gmx_mtop_t *mtop, rvec v[]);

void hl_gaussian(rvec v[], int natoms, real sd, unsigned int seed); /* added by Huiyu Li, 2023 02 07, generate Gaussian */

/*! \brief
 * Remove the center of mass motion in a set of coordinates.
 *
 * \param[out] log  File for printing debug information
 * \param[in]  natoms Number of atoms
 * \param[in]  mass   Atomic masses
 * \param[in]  x      Coordinates
 * \param[out] v      Velocities
 */
void stop_cm(FILE *log, int natoms, real mass[], rvec x[], rvec v[]);

#endif
