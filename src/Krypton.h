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

#ifndef KRYPTON_H
#define KRYPTON_H
#include "Atom.h"

class Krypton : public Atom
{
private:
    static const int Z = 36;
    static const double mass = 84.0;
    
    static const double ljDimer = 7.58; // DOI:10.1080/00268978900101821

    static const double ionization_potentials[Z];

    static const int norbitals = 8;
    static const int ground_state_configuration[norbitals];

public:

    double Eip ( int q ) {
        return ionization_potentials[q];
    }

    double getMass()
    {
        return mass * 1822.88848426455;
    }
    
    double getLJDimer()
    {
        return ljDimer;
    }
};

const int Krypton::ground_state_configuration[Krypton::norbitals] = {2, 2, 6, 2, 6, 10, 2, 6};

const double   Krypton::ionization_potentials[Krypton::Z] = {
    0.514505696435134,0.895240720323411,1.35787210584344,1.9311098860713,
    2.37675119441382,2.88333333333333,4.07930907754502,4.62322675486953,
    8.48393972804116,9.85740536567438,11.3124219037119,12.8740536567438,
    14.3595369349504,16.4163175303197,18.0922454979787,19.8824329290702,
    21.7487688349871,23.5389562660786,28.871407570746,30.6235207644248,
    32.4898566703418
};

#endif // KRYPTON_H
