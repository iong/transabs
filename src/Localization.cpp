/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/
#include <cmath>

#include "Common.h"
#include "Localization.h"

static double angularVelocity (int ei, int aj)
{
    double dx, dy, dz, dvx, dvy, dvz, Lx, Ly, Lz;
    double angvel;

    /** \f[ d\Theta = \frac{\vec{v}\cross\vec{r}}{r^2} dt \f] */
    dvx = vx[ei] - vx[aj];
    dvy = vy[ei] - vy[aj];
    dvz = vz[ei] - vz[aj];

    dx = x[ei] - x[aj];
    dy = y[ei] - y[aj];
    dz = z[ei] - z[aj];

    Lx = dy * dvz - dz * dvy;
    Ly = dz * dvx - dx * dvz;
    Lz = dx * dvy - dy * dvx;

    angvel = (Lx * Lx + Ly * Ly + Lz * Lz) / (dx * dx + dy * dy + dz * dz);

    return angvel;
}


void Localization (double dt)
{
    size_t i;

    for (i = 0; i < Nparticles - Natom; i++)
    {
        if (new_next_atom[i] != next_atom[i])
        {
            if (next_atom[i] >= 0 && revangle[i] > LocalizationAngle) {
                assert(nloc[next_atom[i]]>0);
                nloc[next_atom[i]]--;
            }
            valence[i] = 0;

            next_atom[i] = new_next_atom[i];
            next_atom_dist[i] = new_next_atom_dist[i];
            revangle[i] = 0.0;
        }
        else
        {
            double oldangle = revangle[i];
            revangle[i] += angularVelocity (i+Natom, next_atom[i]) * dt;
            next_atom_dist[i] = sqrt(new_next_atom_dist[i]);

            if (revangle[i] >= LocalizationAngle && oldangle < LocalizationAngle)
            {
                nloc[next_atom[i]]++;
            }
        }
    }
}