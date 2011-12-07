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

#ifndef CUBEHISTOGRAM_H
#define CUBEHISTOGRAM_H

#include <cstddef>
#include <cmath>
#include <armadillo>

using namespace arma;
using namespace std;

class CubeHistogram
{
public:
    u32  nxbins, nybins, nzbins;
    double  xmin, xmax, dx;
    double  ymin, ymax, dy;
    double  zmin, zmax, dz;
    vec     xrange, yrange, zrange;
    cube    bins;

    CubeHistogram () :
            nxbins(0),
            nybins(0),
            nzbins(0),
            xmin (0),
            xmax (0),
            dx (1.0),
            ymin (0),
            ymax (0),
            dy (1.0),
            zmin (0),
            zmax (0),
            dz (1.0)
    {
    }

    CubeHistogram (double _xmin, double _xmax, double _dx) :
            xmin (_xmin),
            xmax (_xmax),
            dx (_dx),
            ymin (_xmin),
            ymax (_xmax),
            dy (_dx),
            zmin (_xmin),
            zmax (_xmax),
            dz (_dx)
    {
        nxbins = round ( (xmax - xmin) / dx);
        nybins = nxbins;
        nzbins = nxbins;

        bins.zeros(nxbins, nybins, nzbins);
    }

    CubeHistogram (double _xmin, double _dx, int _nxbins) :
            nxbins(_nxbins), xmin (_xmin), dx(_dx) 
    {

        xmax = xmin + dx * nxbins;

        ymin = xmin;
        ymax = xmax;
        dy = dx;

        zmin = xmin;
        zmax = xmax;
        dz = dx;

        nybins = nxbins;
        nzbins = nxbins;

        bins.zeros(nxbins, nybins, nzbins);
    }

    CubeHistogram (double _xmin, double _dx, u32 _nxbins,
                   double _ymin, double _dy, u32 _nybins,
                   double _zmin, double _dz, u32 _nzbins) :
            nxbins(_nxbins), nybins(_nybins), nzbins(_nzbins),
            xmin (_xmin), dx(_dx), ymin (_ymin), dy(_dy),
            zmin (_zmin), dz(_dz)
    {
        xmax = xmin + dx * nxbins;
        ymax = ymin + dy * nybins;
        zmax = zmin + dz * nzbins;

        bins.zeros(nxbins, nybins, nzbins);
    }

    CubeHistogram (double _xmin, double _xmax, double _dx,
                   double _ymin, double _ymax, double _dy,
                   double _zmin, double _zmax, double _dz) :
            xmin (_xmin),
            xmax (_xmax),
            dx (_dx),
            ymin (_ymin),
            ymax (_ymax),
            dy (_dy),
            zmin (_zmin),
            zmax (_zmax),
            dz (_dz)
    {
        nxbins =  round ( (xmax - xmin) / dx);
        nybins =  round ( (ymax - ymin) / dy);
        nzbins =  round ( (zmax - zmin) / dz);

        bins.zeros(nxbins, nybins, nzbins);
    }

    bool empty() {
        return nxbins*nybins*nzbins == 0;
    }

    vec get_xrange()
    {
        vec xrange(nxbins);

        for (u32 i = 0; i<nxbins; i++) {
            xrange[i] = xmin + dx*(double)i;
        }

        return xrange;
    }

    void increment (double x, double y, double z)
    {
        u32 i, j, k;

        if (x < xmin || y < ymin || z < zmin
                || x >= xmax || y >= ymax || z >= zmax) return;

        i = floor ( (x - xmin) / dx);
        j = floor ( (y - ymin) / dy);
        k = floor ( (z - zmin) / dz);

        bins(i, j, k) += 1.0;
    }

    void increment (vec &x, vec &y, vec &z)
    {
        vec xi(x.n_rows), yi(y.n_rows), zi(y.n_rows);

        xi = floor( (x-xmin) / dx);
        yi = floor( (y-ymin) / dy);
        zi = floor( (z-zmin) / dz);

        for (u32 k=0; k<x.n_elem; k++) {
            if (xi(k) < 0 || yi(k) < 0 || zi(k) < 0
                    || xi(k) >= nxbins || yi(k) >= nybins || zi(k) >= nzbins)
                continue;

            bins( xi(k), yi(k), zi(k) ) += 1.0;
        }
    }

    template<typename eT>
    void increment (const subview_col<eT> &x, const subview_col<eT> &y, const subview_col<eT> &z)
    {
        const u32 n_elem     = x.n_elem;
        Col<eT> xi(n_elem), yi(n_elem), zi(n_elem);

        xi = floor( (x-xmin) / dx);
        yi = floor( (y-ymin) / dy);
        zi = floor( (z-zmin) / dz);

        for (u32 k=0; k<n_elem; k++) {
            if (xi(k) < 0 || yi(k) < 0 || zi(k) < 0
                    || xi(k) >= nxbins || yi(k) >= nybins || zi(k) >= nzbins)
                continue;

            bins( xi(k), yi(k), zi(k) ) += 1.0;
        }
    }

    void increment (vec &x, vec &y, double z)
    {
        vec xi(x.n_rows), yi(y.n_rows);
        double zi;

        xi = floor( (x-xmin) / dx);
        yi = floor( (y-ymin) / dy);
        zi = floor( (z-zmin) / dz);

        for (u32 k=0; k<x.n_elem; k++) {
            if (xi(k) < 0 || yi(k) < 0 || zi < 0
                    || xi(k) >= nxbins || yi(k) >= nybins || zi >= nzbins)
                continue;

            bins( xi(k), yi(k), zi ) += 1.0;
        }
    }

    template<typename eT>
    void increment (subview_col<eT> &x, subview_col<eT> &y, eT z)
    {
        const u32 n_elem     = x.n_elem;
        Col<eT> xi(n_elem);
        Col<eT> yi(n_elem);
        eT zi;

        xi = floor( (x-xmin) / dx);
        yi = floor( (y-ymin) / dy);
        zi = floor( (z-zmin) / dz);


        for (u32 k=0; k<n_elem; k++) {
            if (xi(k) < 0 || yi(k) < 0 || zi < 0
                    || xi(k) >= nxbins || yi(k) >= nybins || zi >= nzbins)
                continue;

            bins( xi(k), yi(k), zi ) += 1.0;
        }
    }
};

#endif // CUBEHISTOGRAM_H
