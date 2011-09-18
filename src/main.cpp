//
//  main.cpp
//  transabs
//
//  Created by Ionuț Georgescu on 8/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <xmmintrin.h>


#include <boost/filesystem/operations.hpp>
#include <boost/assign/list_inserter.hpp>

#include "Common.h"
#include "Argon.h"
#include "Laser.h"
#include "DirectSum.h"
#include "HDF5IO.h"
#include "Electron.h"
#include "Localization.h"
#include "Statistics.h"


using namespace std;

namespace fs = boost::filesystem;

using namespace boost::assign;

po::variables_map vm;

size_t Natom, Nparticles = 0;
vec    q, mass, x, y, z, vx, vy, vz, fx, fy, fz, phi, next_atom_dist,
new_next_atom_dist, revangle, active_electron_delay;
ivec next_atom, new_next_atom, nlocByRevAngle, nlocByEnergy, valence;
span allAtoms, allElectrons, allParticles;
double  Utot, Eoffset, rcluster;

VectorMap<double>   fDumpVectors, fVectors;
VectorMap<int>  iVectors;


cube quasi_free_hist, valence_hist;
static cube quasi_free_rho, valence_rho;
SimpleHistogram2D radialDist;
size_t histogramNo = 0;
double histogram_dt;

vec ne_free, ne_localizedByEnergy, ne_localizedByRevAngle, ne_total, Ekinavg,
Epotavg, EoffsetAvg, statTime;
mat cm, vcm, qavg;

int Nruns = 1;

Argon   sp;

double t, tstart, tstop, dt, dtSnapshot = 10.0;
double soft_core, LocalizationAngle = 3.0;

GaussianLaserPulse pump;


void process_options (int argc, char * argv[])
{
    try
    {

        po::options_description cmdline ("Command line options");
        cmdline.add_options()
        ("help", "get help")
        ("config", po::value<string>(), "configuration file");

        po::options_description generic ("General options");
        generic.add_options()
        ("pump.wave-length", po::value<double>()->required(), "wave length of the pump laser [nm]")
        ("pump.intensity", po::value<double>()->required(), "intensity of the pump laser [W/cm^2]")
        ("pump.fwhm", po::value<double>()->required(), "pulse length of the pump laser [fs]")
        ("coordinates", po::value<string>()->required(), "initial coordinates of the cluster atoms")
        ("eps", po::value<double>()->required(), "soft-core parameter")
        ("tstart", po::value<double>()->required(), "soft-core parameter")
        ("tstop", po::value<double>()->required(), "soft-core parameter")
        ("dt", po::value<double>()->required(), "soft-core parameter")
        ("snapshot.interval", po::value<double>()->required(), "time interval for data dump")
        ("snapshot.file", po::value<string>()->required(), "dump data to this HDF5 file")
        ("localization_angle", po::value<double>()->required(), "localization angle in multiples of pi")
        ("histogram.rmax", po::value<double>()->required(), "maximum radius of the electron histogram")
        ("histogram.dr", po::value<double>()->required(), "histogram bin width")
        ("histogram.dt", po::value<double>()->required(), "time slice \"width\" of the histogram")
        ("nruns", po::value<int>(), "# of trajectories to average over")
        ;

        cmdline.add (generic);

        po::positional_options_description  p;
        p.add ("coordinates", 1);

        po::store (po::command_line_parser (argc, argv).options (cmdline).positional (p).run(), vm);

        if (vm.count ("config"))
        {
            ifstream cfgin (vm["config"].as<string>().c_str());
            po::store (po::parse_config_file (cfgin, generic), vm);
            cfgin.close();
        }

        po::notify (vm);

        if (vm.count ("help"))
        {
            cout << "\n\t\tAsk Ionut!\n\n" << cmdline << endl;
            exit (EXIT_FAILURE);
        }
    }
    catch (po::required_option& ro)
    {
        cerr << "Required option " << ro.get_option_name() << " is missing.\n";
        exit (EXIT_FAILURE);
    }
    catch (exception& e)
    {
        cerr << "error: " << e.what() << "\n";
        exit (EXIT_FAILURE);
    }

    soft_core = vm["eps"].as<double>();

    tstart = vm["tstart"].as<double>();
    tstop = vm["tstop"].as<double>();
    dt = vm["dt"].as<double>();

    // I_au = c eps_0 E_au^2/2
    // w_au = 2*pi*\hbar*c/E_H
    pump = GaussianLaserPulse (vm["pump.intensity"].as<double>() / 3.50944493695457e+16,
                               45.5633525101396 / vm["pump.wave-length"].as<double>(),
                               vm["pump.fwhm"].as<double>() * 41.3413733524035);
    pump.setZero (2.5 * vm["pump.fwhm"].as<double>() * 41.3413733524035);

    dtSnapshot = vm["snapshot.interval"].as<double>();
    LocalizationAngle = vm["localization_angle"].as<double>() * M_PI;

    radialDist = SimpleHistogram2D (0.0,
                                    vm["histogram.rmax"].as<double>(),
                                    vm["histogram.dr"].as<double>(),
                                    0.0,
                                    11.0,
                                    1.0);

    histogram_dt = vm["histogram.dt"].as<double>();

    if (vm.count ("nruns"))
    {
        Nruns = vm["nruns"].as<int>();
    }

}


template<typename T>
void init_fields (VectorMap<T>& f, int N)
{
    typename VectorMap<T>::iterator j;

    for (j = f.begin(); j != f.end(); j++)
    {
        j->second->zeros (N);
    }
}

template<typename T>
void init_fields (VectorMap<T>& f)
{
    typename VectorMap<T>::iterator j;

    for (j = f.begin(); j != f.end(); j++)
    {
        j->second->zeros ();
    }
}


void resize_vectors ()
{
    VectorMap<double>::iterator jf;
    VectorMap<int>::iterator ji;
    int N;

    N = 1.5 * fVectors.begin()->second->n_rows;

    for (jf = fVectors.begin(); jf != fVectors.end(); jf++)
    {
        jf->second->reshape (N, 1);
    }

    for (ji = iVectors.begin(); ji != iVectors.end(); ji++)
    {
        ji->second->reshape (N, 1);
    }
}

static void collect_statistics (size_t  j)
{
    statTime (j) += t;

    vec Ebody (q.n_rows);

    Ebody = kineticEnergy() + q % phi;

    for (size_t k = Natom; k < Nparticles; k++)
    {
        if (Ebody(k) >= 0)
        {
            ne_free(j) += 1.0;
        }
    }

    ne_localizedByEnergy(j) += sum (nlocByEnergy);

    ne_localizedByRevAngle(j) += sum (nlocByRevAngle);

    ne_total (j) += Nparticles - Natom;

    cm (0, j) += mean (x.rows (0, Nparticles - 1));
    cm (1, j) += mean (y.rows (0, Nparticles - 1));
    cm (2, j) += mean (z.rows (0, Nparticles - 1));

    vcm (0, j) += mean (vx.rows (0, Nparticles - 1));
    vcm (1, j) += mean (vy.rows (0, Nparticles - 1));
    vcm (2, j) += mean (vz.rows (0, Nparticles - 1));

    Ekinavg (j) += sum (kineticEnergy());
    Epotavg (j) += Utot;
    EoffsetAvg (j) += Eoffset;

    qavg.col (j) += q.rows (0, Natom - 1);
}



int main (int argc, char * argv[])
{
    mat rt;
    int ndt, ndtSnapshot, ndtStatistics, nStatSlices, NMaxParticles;

    _mm_setcsr (_MM_MASK_MASK &~ (_MM_MASK_OVERFLOW | _MM_MASK_INVALID | _MM_MASK_DIV_ZERO));

    process_options (argc, argv);

    rt.load (vm["coordinates"].as<string>());
    rt *= 3.84 * 1.89;
    Natom = rt.n_rows;
    NMaxParticles = 5 * Natom;

    ndt = (tstop - tstart) / dt;
    ndtSnapshot = dtSnapshot / dt;
    ndtStatistics = histogram_dt / dt;
    nStatSlices = (tstop - tstart) / histogram_dt + 1;



    insert (fDumpVectors) ("q",  &q) ("mass", &mass)
    ("x", &x) ("y", &y) ("z", &z)
    ("vx", &vx) ("vy", &vy) ("vz", &vz)
    ("fx", &fx) ("fy", &fy) ("fzbre", &fz)
    ("phi", &phi);

    fVectors = fDumpVectors;
    insert (fVectors) ("next_atom_dist", &next_atom_dist)
    ("new_next_atom_dist", &new_next_atom_dist) ("revangle", &revangle);

    insert (iVectors) ("next_atom", &next_atom)
    ("new_next_atom", &new_next_atom) ("valence", &valence);

    init_fields (fVectors, NMaxParticles);
    init_fields (iVectors, NMaxParticles);

    allAtoms = span (0, Natom);
    allParticles = span (0, Natom);
    allElectrons = span (0);

    map<string, mat *> statFields;
    insert (statFields) ("ne_free", &ne_free)
    ("ne_localizedByEnergy", &ne_localizedByEnergy)
    ("ne_localizedByRevAngle", &ne_localizedByRevAngle)
    ("ne_total", &ne_total)
    ("Ekinavg", &Ekinavg)
    ("Epotavg", &Epotavg)
    ("EoffsetAvg", &EoffsetAvg)
    ("time", &statTime);

    for (map<string, mat *>::iterator j = statFields.begin();
            j != statFields.end(); j++)
    {
        j->second->zeros (nStatSlices);
    }

    insert (statFields) ("cm", &cm) ("vcm", &vcm) ("qavg", &qavg);

    cm.zeros (3, nStatSlices);
    vcm = cm;
    qavg.zeros (Natom, nStatSlices);

    quasi_free_hist.zeros (radialDist.nxbins, radialDist.nybins, nStatSlices);

    valence_hist = quasi_free_hist;


    HDF5IO h5dump (vm["snapshot.file"].as<string>());

    for (int run = 0; run < Nruns; run++)
    {
        init_fields (fVectors);
        init_fields (iVectors);

        q.zeros();
        mass.ones();
        next_atom.fill (-2);
        nlocByRevAngle.zeros (Natom);
        nlocByEnergy.zeros (Natom);

        x.rows (0, Natom - 1) = rt.col (0);
        y.rows (0, Natom - 1) = rt.col (1);
        z.rows (0, Natom - 1) = rt.col (2);

        Nparticles = Natom;


        mass.rows (0, Natom - 1).fill (sp.getMass());

        for (size_t i = 0; i < Natom; i++)
        {
            create_electron (&sp, i);
        }
        active_electron_delay.set_size(Natom);
        active_electron_delay.fill(-1.0);

        DirectSum();

        int ndtRClusterUpdate = 2.0 * 41 / dt; // every 2fs

        for (int i = 0; i <= ndt; i++, t = tstart + dt * i)
        {
            vec invmass = 1.0 / mass;

            if (i % ndtRClusterUpdate == 0)
            {
                rcluster = clusterRadius();
            }

            if (i % ndtStatistics == 0)
            {
                histogramNo = i / ndtStatistics;
                collect_statistics (histogramNo);
            }

            if (i % ndtSnapshot == 0)
            {
                h5dump.newGroup();
                h5dump.add (fDumpVectors, Nparticles);
                h5dump.add ("nlocByEnergy", nlocByEnergy, Natom);
                h5dump.add ("nlocByRevAngle", nlocByRevAngle, Natom);
                h5dump.add ("next_atom", next_atom, Nparticles);
                h5dump.add ("revangle", revangle, Nparticles);
                h5dump.closeGroup();
            }

            vx += (0.5 * dt) * fx % q % invmass;

            vy += (0.5 * dt) * fy % q % invmass;
            vz += (0.5 * dt) * (fz + pump.field (t)) % q % invmass;

            x += dt * vx;
            y += dt * vy;
            z += dt * vz;

            DirectSum();

            vx += (0.5 * dt) * fx % q % invmass;
            vy += (0.5 * dt) * fy % q % invmass;
            vz += (0.5 * dt) * (fz + pump.field (t)) % q % invmass;

            Localization (dt);

            for (size_t j = 0; j < Natom; j++)
            {
                if (nlocByEnergy(j) > 0) continue;

                create_electron (&sp, j);
/*
                if (active_electron_delay[j] == -1.0)
                {
		    active_electron_delay[j] = 20; //2.0 * M_PI * q[j] / pow (sp.Eip (q[j] - 1), 1.5);
                }
                else if (active_electron_delay[j] > 0.0)
                {
		    active_electron_delay[j] -= dt;
                }
                else
                {
		    create_electron (&sp, j);
		    active_electron_delay[j] = -1.0;
                }
*/
	    }
        }
    }

    quasi_free_hist *=  dt / (histogram_dt * (double) Natom * Nruns);

    valence_hist *=  dt / (histogram_dt * (double) Natom * Nruns);

    vec shellVolumes, rbins (radialDist.get_xrange()), rmax;
    rmax = radialDist.xmax;

    shellVolumes = pow (join_cols (rbins.rows (1, rbins.n_elem - 1), rmax), 3) -
                   pow (rbins, 3);
    shellVolumes *= 4.0 * M_PI / 3.0;


    quasi_free_rho = quasi_free_hist;
    valence_rho = valence_hist;

    for (size_t j = 0; j < quasi_free_rho.n_slices; j++)
    {
        for (size_t k = 0; k < quasi_free_rho.n_cols; k++)
        {
            quasi_free_rho.slice (j).col (k) /= shellVolumes;
            valence_rho.slice (j).col (k) /= shellVolumes;
        }
    }

    vec laserField (statTime.n_rows);

    for (size_t j = 0; j < statTime.n_rows; j++)
    {
        laserField (j) = pump.field (statTime (j));
    }

    insert (statFields) ("laserField", &laserField);


    h5dump.addStatistics ("quasi_free_hist", quasi_free_hist);

    h5dump.addStatistics ("quasi_free_rho", quasi_free_rho);
    h5dump.addStatistics ("valence_hist", valence_hist);
    h5dump.addStatistics ("valence_rho", valence_rho);
    h5dump.addStatistics ("rbins", rbins);

    for (map<string, mat *>::iterator j = statFields.begin(); j != statFields.end(); j++)
    {
        * (j->second) /= (double) Nruns;
        h5dump.addStatistics (j->first.c_str(), * (j->second));
    }


    return 0;
}


