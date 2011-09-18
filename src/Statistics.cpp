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
#include "Statistics.h"

void atomsNextAtoms()
{
    vec	rsq(Natom);
    
    next_atom_dist.rows(0, Natom -1).fill(INFINITY);
    
    
    for (size_t i = Natom - 1; i > 0; i--) {
	span ispan(0, i - 1);
	
	rsq.rows(0, i-1) = square(x.rows(0, i-1) - x(i)) + square(y.rows(0, i-1) - y(i)) + square(z.rows(0, i-1) - z(i));

	for (size_t j = 0; j<i; j++) {
	    if (rsq(j) < next_atom_dist(j)) {
		next_atom(j) = i;
		next_atom_dist(j) = rsq(j);
	    }
	    if (rsq(j) < next_atom_dist(i)) {
		next_atom(i) = j;
		next_atom_dist(i) = rsq(j);
	    }
	}
    }
    next_atom_dist.rows(0, Natom -1) = sqrt(next_atom_dist.rows(0, Natom -1));
}

double clusterRadius()
{
    vec r(Natom);
    double	rcluster;
    
    atomsNextAtoms();
    
    r = sqrt(square(x.rows(0, Natom -1)) + square(y.rows(0, Natom -1)) + square(z.rows(0, Natom -1)));
    r += next_atom_dist.rows(0, Natom - 1);
    
    rcluster = max(r);
    
    return rcluster;
}

vec kineticEnergy()
{
    vec	Ekin;
    Ekin =  0.5 * (square (vx) + square (vy) + square (vz)) % mass;
    return Ekin;
}


double kineticEnergy(size_t k)
{
    return 0.5 * mass(k) * (vx(k) * vx(k) + vy(k) * vy(k) + vz(k) * vz(k));
}