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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <boost/assign/list_inserter.hpp>

#include "Statistics.h"

#include "Common.h"
#include "CubeHistogram.h"
#include "HDF5IO.h"
#include "Laser.h"
#include "VectorMap.h"

using namespace boost::assign;

extern LaserPulse *pump;

int statSlabLength, statNSlabs;

CubeHistogram   electronDensity;

CubeHistogram hist_qfe, hist_valence, hist_localized, hist_delocalized;

cube rho_localized, rho_delocalized;
Cube<int>    radDistNorm;

map<string, mat *> statFields;

vec ne_free, ne_localizedByEnergy, ne_localizedByRevAngle, ne_total, Ekinavg,
    Epotavg, EoffsetAvg, statTime;
mat cm, vcm, qavg, rcmQuasiFreeElectrons, vcmQuasiFreeElectrons,
    rcmClusterElectrons, vcmClusterElectrons;

static int qMax;
static double histogram_dt;

void initStats()
{

    histogram_dt = vm["histogram.dt"].as<double>();

    statSlabLength = histogram_dt / dt;
    statNSlabs = (tstop - tstart) / histogram_dt + 1;

    insert (statFields) ("ne_free", &ne_free)
    ("ne_localizedByEnergy", &ne_localizedByEnergy)
    ("ne_localizedByRevAngle", &ne_localizedByRevAngle)
    ("ne_total", &ne_total)
    ("Ekinavg", &Ekinavg)
    ("Epotavg", &Epotavg)
    ("EoffsetAvg", &EoffsetAvg)
    ("time", &statTime);

    for (map<string, mat *>::iterator j = statFields.begin(); j != statFields.end(); j++)
    {
        j->second->zeros (statNSlabs);
    }

    insert (statFields) ("cm", &cm) ("vcm", &vcm) ("qavg", &qavg)
    ("rcmQuasiFreeElectrons", &rcmQuasiFreeElectrons)
    ("vcmQuasiFreeElectrons", &vcmQuasiFreeElectrons)
    ("rcmClusterElectrons", &rcmClusterElectrons)
    ("vcmClusterElectrons", &vcmClusterElectrons);

    cm.zeros (3, statNSlabs);
    vcm = cm;
    rcmQuasiFreeElectrons = cm;
    vcmQuasiFreeElectrons = cm;
    rcmClusterElectrons = cm;
    vcmClusterElectrons = cm;

    qavg.zeros (Natom, statNSlabs);

    qMax = vm["histogram.maximumCharge"].as<int>();

    hist_localized = CubeHistogram(0.0, vm["histogram.rmax"].as<double>(), vm["histogram.dr"].as<double>(),
                             0, Natom*(qMax   + 1), 1.0,
                             tstart, tstop, histogram_dt);
    hist_delocalized = hist_localized;

    radDistNorm.zeros(Natom, qMax+1, hist_localized.nzbins);

    if (vm.count ("grid.dx"))
    {
        double dx = vm["grid.dx"].as<double>();
        int gridsize = vm["grid.size"].as<int>();

        electronDensity  = CubeHistogram (-0.5*dx*gridsize, dx, gridsize,
                                          -0.5*dx*gridsize, dx, gridsize,
                                          -0.5*dx*gridsize, dx, gridsize*statNSlabs);
    }

}

void atomsNextAtoms()
{
    vec rsq(Natom);

    next_atom_dist.rows(0, Natom -1).fill(INFINITY);


    for (size_t i = Natom - 1; i > 0; i--) {
        rsq.rows(0, i-1) = square(x.rows(0, i-1) - x(i))
                           + square(y.rows(0, i-1) - y(i))
                           + square(z.rows(0, i-1) - z(i));

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
    double  rcluster;

    atomsNextAtoms();

    r = sqrt(square(x.rows(0, Natom -1)) + square(y.rows(0, Natom -1)) + square(z.rows(0, Natom -1)));
    r += next_atom_dist.rows(0, Natom - 1);

    rcluster = max(r);

    return rcluster;
}

vec kineticEnergy()
{
    vec Ekin;
    Ekin =  0.5 * (square (vx) + square (vy) + square (vz)) % mass;
    return Ekin;
}


double kineticEnergy(size_t k)
{
    return 0.5 * mass(k) * (vx(k) * vx(k) + vy(k) * vy(k) + vz(k) * vz(k));
}

static void sliceStats (size_t  slabNo)
{
    int nQuasiFreeElectrons, nClusterElectrons;
    statTime (slabNo) += t;

    vec Ebody (q.n_rows), rsq(q.n_rows);

    Ebody = kineticEnergy() + q % phi;
    rsq = square(x) + square(y) + square(z);

    double clusterLimitSq = 1.25*1.25*rcluster*rcluster;

    vec::fixed<3> rcmQuasiFreeElectrons_, vcmQuasiFreeElectrons_,
    rcmClusterElectrons_, vcmClusterElectrons_;

    rcmQuasiFreeElectrons_.zeros();
    vcmQuasiFreeElectrons_.zeros();
    rcmClusterElectrons_.zeros();
    vcmClusterElectrons_.zeros();
    nQuasiFreeElectrons = 0;
    nClusterElectrons = 0;

    for (size_t k = Natom; k < Nparticles; k++)
    {
        if (Ebody(k) >= 0)
        {
            ne_free(slabNo) += 1.0;
        }
        else {
            rcmQuasiFreeElectrons_[0] += x[k];
            rcmQuasiFreeElectrons_[1] += y[k];
            rcmQuasiFreeElectrons_[2] += z[k];

            vcmQuasiFreeElectrons_[0] += vx[k];
            vcmQuasiFreeElectrons_[1] += vy[k];
            vcmQuasiFreeElectrons_[2] += vz[k];

            nQuasiFreeElectrons++;
        }

        if (rsq [k] <= clusterLimitSq) {
            rcmClusterElectrons_[0] += x[k];
            rcmClusterElectrons_[1] += y[k];
            rcmClusterElectrons_[2] += z[k];

            vcmClusterElectrons_[0] += vx[k];
            vcmClusterElectrons_[1] += vy[k];
            vcmClusterElectrons_[2] += vz[k];

            nClusterElectrons++;
        }
    }

    ne_localizedByEnergy(slabNo) += sum (nlocByEnergy); // A lot of double counting?

    ne_localizedByRevAngle(slabNo) += sum (nlocByRevAngle);

    ne_total (slabNo) += Nparticles - Natom;

    cm (0, slabNo) += mean (x.rows (0, Nparticles - 1));
    cm (1, slabNo) += mean (y.rows (0, Nparticles - 1));
    cm (2, slabNo) += mean (z.rows (0, Nparticles - 1));

    vcm (0, slabNo) += mean (vx.rows (0, Nparticles - 1));
    vcm (1, slabNo) += mean (vy.rows (0, Nparticles - 1));
    vcm (2, slabNo) += mean (vz.rows (0, Nparticles - 1));

    rcmQuasiFreeElectrons.col(slabNo) += rcmQuasiFreeElectrons_/nQuasiFreeElectrons;
    vcmQuasiFreeElectrons.col(slabNo) += vcmQuasiFreeElectrons_/nQuasiFreeElectrons;

    rcmClusterElectrons.col(slabNo) += rcmClusterElectrons_/nClusterElectrons;
    vcmClusterElectrons.col(slabNo) += vcmClusterElectrons_/nClusterElectrons;

    Ekinavg (slabNo) += sum (kineticEnergy());
    Epotavg (slabNo) += Utot;
    EoffsetAvg (slabNo) += Eoffset;

    qavg.col (slabNo) += q.rows (0, Natom - 1);
}

void collectStats(int tslice)
{
    vec shiftedZ(Nparticles - Natom);

    int slabNo = tslice / statSlabLength;

    if (!electronDensity.empty()) {
        shiftedZ = z.rows(Natom, Nparticles-1)
                   + slabNo * (electronDensity.nzbins/statNSlabs)
                   * electronDensity.dz;

        electronDensity.increment(x.rows(Natom, Nparticles-1),
                                  y.rows(Natom, Nparticles-1),
                                  shiftedZ.subvec(span::all));
    }

    if (tslice % statSlabLength != 0) {
        return;
    }

    sliceStats(tslice / statSlabLength);
}

void incrementRadialDistributions(size_t ie, vec &atomDist, vec &aid)
{

    if (revangle[ie] >= LocalizationAngle)
    {
        hist_localized.increment(atomDist, aid, t);
    }
    else
    {
        hist_delocalized.increment(atomDist, aid, t);
    }
}

void incrementRadialDistributionsQ(size_t ie, vec &atomDist, vec &realCharge)
{
    vec idx(Natom);

    int tslice = floor(t - tstart) / histogram_dt;
    for (int i=0; i++; i < realCharge.n_elem) {
        idx(i) = min(realCharge(i), (double)qMax) * Natom + i;
        radDistNorm(i, (int)realCharge(i), tslice) += 1;
    }

    if (revangle[ie] >= LocalizationAngle)
    {
        hist_localized.increment(atomDist, idx, t);
    }
    else
    {
        hist_delocalized.increment(atomDist, idx, t);
    }
}

#ifdef HAVE_MPI
template <typename ms>
static void mpiReduce(ms &c, MPI_Datatype dtype)
{
    if (mpiRank != 0) {
        MPI_Reduce(c.memptr(), NULL, c.n_elem, dtype, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Reduce(MPI_IN_PLACE, c.memptr(), c.n_elem, dtype, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}

void centralizeStats()
{
    for (map<string, mat *>::iterator j = statFields.begin(); j != statFields.end(); j++)
    {
        mpiReduce<mat>(* (j->second), MPI_DOUBLE);
    }

    mpiReduce<cube>(hist_localized.bins, MPI_DOUBLE);
    mpiReduce<cube>(hist_delocalized.bins, MPI_DOUBLE);
    mpiReduce<cube>(electronDensity.bins, MPI_DOUBLE);
    mpiReduce<Cube <int> >(radDistNorm, MPI_INT);
}
#endif

void normalizeStats(int nRuns)
{
    // Normalize all histograms
    for (map<string, mat *>::iterator j = statFields.begin(); j != statFields.end(); j++)
    {
        * (j->second) /= (double) nRuns;
    }

    double norm;
    for (size_t j = 0; j < radDistNorm.n_slices; j++) {
        for (size_t k = 0; k < radDistNorm.n_cols; k++) {
            for (size_t l = 0; l < radDistNorm.n_cols; l++) {
		norm = radDistNorm(l, k, j);
                if ( radDistNorm(l, k, j) == 0.0) continue;

                hist_localized.bins.slice(j).col(k*Natom + l) /= norm;
                hist_delocalized.bins.slice(j).col(k*Natom + l) /= norm;
            }
        }
    }

    if (!electronDensity.empty()) {
        electronDensity.bins *= 1.0 / (double)(statSlabLength * nRuns);
    }
    // Convert histograms to densities
    vec shellVolumes, rbins (hist_localized.get_xrange()), rmax;
    rmax = hist_localized.xmax;

    shellVolumes = join_cols (rbins.rows (1, rbins.n_elem - 1), rmax);
    shellVolumes = 4.0/3.0 * M_PI * ( pow(shellVolumes, 3) - pow(rbins, 3) );

    rho_localized = hist_localized.bins;
    rho_delocalized = hist_delocalized.bins;

    for (size_t j = 0; j < rho_localized.n_slices; j++)
    {
        for (size_t k = 0; k < rho_localized.n_cols; k++)
        {
            rho_localized.slice (j).col (k) /= shellVolumes;
            rho_delocalized.slice (j).col (k) /= shellVolumes;
        }
    }
}



void dumpStats(HDF5IO& h5dump)
{
    // laser field for comparison purposes
    vec laserField (statTime.n_rows);

    for (size_t j = 0; j < statTime.n_rows; j++)
    {
        laserField (j) = pump->field (statTime (j));
    }

    insert (statFields) ("laserField", &laserField);


    vec rbins(hist_qfe.get_xrange());
    h5dump.addCubeSeries1 ("rho_localized", rho_localized, qMax + 1);
    h5dump.addCubeSeries1 ("rho_delocalized", rho_delocalized, qMax + 1);
    h5dump.addField ("radDistNorm", radDistNorm);
    h5dump.addField ("rbins", rbins);
    if (!electronDensity.empty()) {
        h5dump.addCubeSeries2("electronDensity", electronDensity.bins,
                                  statNSlabs);
    }

    for (map<string, mat *>::iterator j = statFields.begin(); j != statFields.end(); j++)
    {
        h5dump.addField (j->first.c_str(), * (j->second));
    }
}
