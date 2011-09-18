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

#ifndef SIMPLEHISTOGRAM2D_H
#define SIMPLEHISTOGRAM2D_H

#include <cmath>
#include <armadillo>

using namespace arma;
using namespace std;

class SimpleHistogram2D
{
    public:
    size_t nxbins, nybins;
    double xmin, xmax, dx;
    double ymin, ymax, dy;
    vec xrange, yrange;
    
    SimpleHistogram2D () :
            xmin (0),
            xmax (0),
            dx (1.0),
            nxbins(0),
            ymin (0),
            ymax (0),
            dy (1.0),
            nybins(0)
    {
    }

    SimpleHistogram2D (double _xmin, double _xmax, double _dx,
		     double _ymin, double _ymax, double _dy) :
            xmin (_xmin),
            xmax (_xmax),
            dx (_dx),
            ymin (_ymin),
            ymax (_ymax),
            dy (_dy)
    {
	nxbins = (size_t) round ( (xmax - xmin) / dx);
	nybins = (size_t) round ( (ymax - ymin) / dy);
    }
    
    vec get_xrange()
    {
	vec xrange(nxbins);

	for (size_t i = 0; i<nxbins; i++) {
	    xrange[i] = xmin + dx*(double)i;
	}
	
	return xrange;
    }
    
    vec get_yrange()
    {
	vec yrange(nybins);

	for (size_t i = 0; i<nybins; i++) {
	    yrange[i] = ymin + dy*(double)i;
	}
	
	return yrange;
    }

    void increment (mat &m, double x, double y)
    {
        u32 i, j;

	if (x < xmin || y < ymin || x >= xmax || y >= ymax) return;

        i = floor ( (x - xmin) / dx);
	j = floor ( (y - ymin) / dy);

        m(i, j) += 1.0;
    }
    
    void increment (mat &m, vec &x, vec &y)
    {
        vec xi(x.n_rows), yi(y.n_rows);

	xi = floor( (x-xmin) / dx);
	yi = floor( (y-ymin) / dy);
	
	for (size_t k=0; k<x.n_elem; k++) {
	    if (xi(k) < 0 || yi(k) < 0 || xi(k) >= nxbins || yi(k) >= nybins) continue;

	    m( xi(k), yi(k) ) += 1.0;
	}
    }
};

#endif // SIMPLEHISTOGRAM2D_H
