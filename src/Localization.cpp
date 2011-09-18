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
#include "Statistics.h"

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
    double Erel;
    size_t	i;
    
    nlocByEnergy.zeros();

    for (i = Natom; i < Nparticles; i++)
    {
	if (new_next_atom(i) > 0 && new_next_atom_dist(i) > next_atom_dist(new_next_atom(i))) {
	    new_next_atom(i) = -1;
	}
	
	
	
        if (new_next_atom(i) != next_atom(i))
        {
            if (next_atom(i) >= 0 && revangle(i) > LocalizationAngle) {
                nlocByRevAngle(next_atom(i))--;
            }
            valence(i) = 0;
            revangle(i) = 0.0;
	    
	    next_atom(i) = new_next_atom(i);
        }
        else if (next_atom(i) >= 0)
        {
            double oldangle = revangle(i);
            revangle(i) += angularVelocity (i, next_atom(i)) * dt;

            if (revangle(i) >= LocalizationAngle && oldangle < LocalizationAngle)
            {
                nlocByRevAngle(next_atom(i))++;
            }
        }
        
        next_atom_dist(i) = new_next_atom_dist(i);
	next_atom(i) = new_next_atom(i);
	// nothing to do: new_next_atom = next_atom = -1
        if (next_atom(i) < 0) {
	    continue;
	}
	
	Erel = kineticEnergy(i) + q(i)*q(next_atom(i))
	    / sqrt(next_atom_dist(i)*next_atom_dist(i) + soft_core*soft_core);
	
	if (Erel < 0) {
	    nlocByEnergy(next_atom(i))++;
	}
    }
}

void LocalizationEnergy()
{
    double Erel;
    size_t	i;
    
    nlocByEnergy.zeros();
    
    for (i = Natom; i < Nparticles; i++)
    {
	next_atom(i) = new_next_atom(i);
	next_atom_dist(i) = new_next_atom_dist(i);
	
	if (next_atom(i) == -1) {
	    continue;
	}
	else if (next_atom(i) > 0 && next_atom_dist(i) > next_atom_dist(next_atom(i))) {
	    next_atom(i) = -1;
	    continue;
	}
        
        Erel = kineticEnergy(i) + q(i)*q(next_atom(i))
	    / sqrt(next_atom_dist(i)*next_atom_dist(i) + soft_core*soft_core);
	
	if (Erel < 0) {
	    nlocByEnergy(next_atom(i))++;
	}
    }
}