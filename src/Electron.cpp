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

#include "Common.h"
#include "Electron.h"
#include "Utils.h"

void create_electron(Atom *a, int aid)
{
    int	eid;

    if ((q[aid] + 1) / soft_core < a->Eip(q[aid])) {
	cerr << "Can't ionize q["<<aid<<"] = "<<q[aid]<<"; (q+1)/a = ";
	cerr << (q[aid] + 1) / soft_core << " < ";
	cerr << " Eip = "<< a->Eip(q[aid]) << endl;
	return;
    }
    
    if (Nparticles + 1 > q.n_elem ){
	resize_vectors();
    }
    
    eid = Nparticles++;

    q[eid] = -1.0;
    mass[eid] = 1.0;

    x[eid] = x[aid];
    y[eid] = y[aid];
    z[eid] = z[aid];
    
    double pol[3];
    pol [0] = sqrt(2.0 * ( (q[aid] + 1) / soft_core - a->Eip(q[aid]) ));
    pol [1] = acos(2.0 * myrand() - 1.0);
    pol [2] = 2.0 * M_PI * myrand();
    polar2cart(pol, vx[eid], vy[eid], vz[eid]);
    
    Eoffset += -a->Eip(q[aid]);
    
    valence[eid] = 1;
    next_atom[eid] = aid;
    next_atom_dist[eid] = 0.0;
    revangle[eid] = LocalizationAngle;
    
    
    q[aid] += 1.0;
    if (q[aid] > 10) {
        cerr << "Stop\n";
    }
    nlocByEnergy[aid]++;
    nlocByRevAngle[aid]++;
}


