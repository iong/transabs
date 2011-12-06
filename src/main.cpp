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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <boost/filesystem/operations.hpp>
#include <boost/assign/list_inserter.hpp>

#include "Common.h"
#include "Krypton.h"
#include "Laser.h"
#include "DirectSum.h"
#include "HDF5IO.h"
#include "Electron.h"
#include "Localization.h"
#include "Statistics.h"
#include "CubeHistogram.h"
#include "VectorMap.h"
#include "Utils.h"


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

int mpiNProcs = 1;
int mpiRank = 0;
int Nruns = 1;

Krypton   sp;

double t, tstart, tstop, dt, dtSnapshot = 10.0;
double soft_core, LocalizationAngle = 3.0;

LaserPulse *pump;

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
        ("pump.phase", po::value<double>(), "pulse phase of the pump laser [multiples of pi]")
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
        ("grid.dx", po::value<double>(), "Electron density grid : dx")
        ("grid.size", po::value<int>(), "Electron density grid : size")
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
    pump = new GaussianLaserPulse (vm["pump.intensity"].as<double>() / 3.50944493695457e+16,
                                   45.5633525101396 / vm["pump.wave-length"].as<double>(),
                                   vm["pump.fwhm"].as<double>() * 41.3413733524035);
    pump->setZero (2.5 * vm["pump.fwhm"].as<double>() * 41.3413733524035);
    if (vm.count ("pump.phase")) {
        pump->setPhase(vm["pump.phase"].as<double>());
    }

    dtSnapshot = vm["snapshot.interval"].as<double>();
    LocalizationAngle = vm["localization_angle"].as<double>() * M_PI;

    if (vm.count ("nruns"))
    {
        Nruns = vm["nruns"].as<int>();
    }

}


void initParticleFields (int NMaxParticles)
{
    typename VectorMap<int>::iterator       ji;
    typename VectorMap<double>::iterator    jf;

    if (fVectors.empty()) {
        insert (fDumpVectors) ("q",  &q) ("mass", &mass)
        ("x", &x) ("y", &y) ("z", &z)
        ("vx", &vx) ("vy", &vy) ("vz", &vz)
        ("fx", &fx) ("fy", &fy) ("fz", &fz)
        ("phi", &phi);

        fVectors = fDumpVectors;
        insert (fVectors) ("next_atom_dist", &next_atom_dist)
        ("new_next_atom_dist", &new_next_atom_dist)
        ("revangle", &revangle);

        insert (iVectors) ("next_atom", &next_atom)
        ("new_next_atom", &new_next_atom)
        ("valence", &valence);
    }

    for (jf = fVectors.begin(); jf != fVectors.end(); jf++)
    {
        jf->second->zeros (NMaxParticles);
    }
    for (ji = iVectors.begin(); ji != iVectors.end(); ji++)
    {
        ji->second->zeros (NMaxParticles);
    }

    allAtoms = span (0, Natom);
    allParticles = span (0, Natom);
    allElectrons = span (0);
}


void resize_vectors ()
{
    VectorMap<double>::iterator jf;
    VectorMap<int>::iterator ji;
    int N, oldN;

    oldN = fVectors.begin()->second->n_rows;
    N = 1.5 * oldN;

    for (jf = fVectors.begin(); jf != fVectors.end(); jf++)
    {
        jf->second->reshape (N, 1);
    }

    for (ji = iVectors.begin(); ji != iVectors.end(); ji++)
    {
        ji->second->reshape (N, 1);
    }
    
    mass.rows(oldN, N-1).fill(1.0);
}



int main (int argc, char * argv[])
{
    mat rt;
    int ndt, ndtSnapshot, NMaxParticles;

    _mm_setcsr (_MM_MASK_MASK &~ (_MM_MASK_OVERFLOW | _MM_MASK_INVALID | _MM_MASK_DIV_ZERO));

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiNProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif

    process_options (argc, argv);

    if (mpiNProcs > Nruns) {
        cerr << "Too many processes("<<mpiNProcs<<")! Only " << Nruns;
        cerr << " needed!\n";
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    int localNruns = (Nruns + mpiNProcs - 1) / mpiNProcs;
    if (mpiRank == mpiNProcs - 1) {
        localNruns = Nruns - (mpiNProcs - 1) * localNruns;
    }

    rt.load (vm["coordinates"].as<string>());
    rt *= sp.getLJDimer();;
    Natom = rt.n_rows;
    NMaxParticles = 5 * Natom;

    initStats();
    init_rng(mpiRank*1000000);

    ndt = (tstop - tstart) / dt;
    ndtSnapshot = dtSnapshot / dt;

    HDF5IO h5dump;
    if (mpiRank == 0) {
        h5dump.open(vm["snapshot.file"].as<string>());
        cout << "HDF5 fid: "<<h5dump.getFID() << endl;
    }
    for (int run = 0; run < localNruns; run++)
    {
        NMaxParticles = NMaxParticles > q.n_elem ? NMaxParticles : q.n_elem;
        initParticleFields(NMaxParticles);

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

            collectStats (i);

            if (mpiRank == 0 && i % ndtSnapshot == 0)
            {
                h5dump.newSnapshot();
                h5dump.addSnapshotFields (fDumpVectors, Nparticles);
                h5dump.addSnapshotField ("nlocByEnergy", nlocByEnergy, Natom);
                h5dump.addSnapshotField ("nlocByRevAngle", nlocByRevAngle, Natom);
                h5dump.addSnapshotField ("next_atom", next_atom, Nparticles);
                h5dump.addSnapshotField ("revangle", revangle, Nparticles);
                h5dump.closeSnapshot();
            }

            vx += (0.5 * dt) * fx % q % invmass;
            vy += (0.5 * dt) * fy % q % invmass;
            vz += (0.5 * dt) * (fz + pump->field (t)) % q % invmass;

            x += dt * vx;
            y += dt * vy;
            z += dt * vz;

            DirectSum();

            vx += (0.5 * dt) * fx % q % invmass;
            vy += (0.5 * dt) * fy % q % invmass;
            vz += (0.5 * dt) * (fz + pump->field (t)) % q % invmass;

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
    
    cout << "ll: " << mpiRank << endl;

#ifdef HAVE_MPI
    centralizeStats();
#endif

    if (mpiRank == 0) {
        normalizeStats(Nruns);
        cout << "HDF5 fid: "<<h5dump.getFID() << endl;
        dumpStats(h5dump);
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}


