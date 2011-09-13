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

using namespace std;

namespace fs = boost::filesystem;

using namespace boost::assign;

po::variables_map vm;

size_t Natom, Nparticles=0;
vec    q, mass, x, y, z, vx, vy, vz, fx, fy, fz, phi, next_atom_dist,
new_next_atom_dist, revangle, next_atom_t0;
ivec next_atom, new_next_atom, nloc, valence;
double  Utot, Eoffset;

VectorMap<double>   fDumpVectors, fVectors;
VectorMap<int>  iVectors;


mat quasi_free_hist, valence_hist;
double histogram_rmax;
double histogram_dr, histogram_dt;
size_t histogramNo=0;
vec histogram_bins, histogram_norm;

vec nfree_electrons, nlocalized_electrons, ntot_electrons, Ekinavg, Epotavg,
    EoffsetAvg, statTime;
mat cm, vcm;

int Nruns=1;

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
    
    histogram_rmax = vm["histogram.rmax"].as<double>();
    histogram_dr = vm["histogram.dr"].as<double>();
    histogram_dt = vm["histogram.dt"].as<double>();
    
    if (vm.count("nruns")) {
        Nruns = vm["nruns"].as<int>();
    }

}


double kinetic_energy()
{
    vec ret (1);

    ret.rows (0, 0) =  0.5 * sum ( (square (vx) + square (vy) + square (vz)) % mass);

    return ret (0);
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
    int	N;
    
    N = 1.5*fVectors.begin()->second->n_rows;

    for (jf = fVectors.begin(); jf != fVectors.end(); jf++)
    {
        jf->second->reshape (N, 1);
    }
    
    for (ji = iVectors.begin(); ji != iVectors.end(); ji++)
    {
        ji->second->reshape (N, 1);
    }
}



int main (int argc, char * argv[])
{
    mat rt;
    int ndt, ndtSnapshot, NMaxParticles;
    
    _mm_setcsr( _MM_MASK_MASK &~ (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );

    process_options (argc, argv);

    rt.load (vm["coordinates"].as<string>());
    rt *= 3.84 * 1.89;
    Natom = rt.n_rows;
    NMaxParticles = 5 * Natom;


    
    insert (fDumpVectors) ("q",  &q) ("mass", &mass)
        ("x", &x) ("y", &y) ("z", &z)
        ("vx", &vx) ("vy", &vy) ("vz", &vz)
        ("fx", &fx) ("fy", &fy) ("fz", &fz)
        ("phi", &phi);

    fVectors = fDumpVectors;
    insert (fVectors) ("next_atom_dist", &next_atom_dist)
        ("new_next_atom_dist", &new_next_atom_dist) ("revangle", &revangle);

    insert (iVectors) ("next_atom", &next_atom)
        ("new_next_atom", &new_next_atom)("valence", &valence);
    
    init_fields (fVectors, NMaxParticles);
    init_fields (iVectors, NMaxParticles);
    
    map<string, mat *> statFields;
    insert(statFields)("nfree_electrons", &nfree_electrons)
            ("nlocalized_electrons", &nlocalized_electrons)
            ("ntot_electrons", &ntot_electrons)
            ("Ekinavg", &Ekinavg)
            ("Epotavg", &Epotavg)
            ("EoffsetAvg", &EoffsetAvg)
            ("time", &statTime);
    
    for (map<string, mat *>::iterator j=statFields.begin();
         j != statFields.end(); j++) {
        
        j->second->zeros((int)floor((tstop-tstart)/histogram_dt) + 1);
    
    }
    
    insert(statFields)("cm", &cm)("vcm", &vcm);
    cm.zeros(3, nfree_electrons.n_rows);
    vcm = cm;
    
    histogram_bins.zeros(floor(histogram_rmax/histogram_dr));
    for (size_t j=0; j<histogram_bins.n_rows; j++) {
        histogram_bins[j] = histogram_dr*j;
    }
    quasi_free_hist.zeros(histogram_bins.n_rows, nfree_electrons.n_rows);
    valence_hist = quasi_free_hist;  
    
    ndt = (tstop - tstart) / dt;
    ndtSnapshot = dtSnapshot / dt;
    HDF5IO h5dump (vm["snapshot.file"].as<string>());
    for (int run=0; run < Nruns; run++) {
        init_fields (fVectors);
        init_fields (iVectors);

        q.zeros();
        mass.ones();
        next_atom.fill (-2);
        nloc.zeros (Natom);

        x.rows (0, Natom - 1) = rt.col (0);
        y.rows (0, Natom - 1) = rt.col (1);
        z.rows (0, Natom - 1) = rt.col (2);
        
        Nparticles = Natom;


        mass.rows (0, Natom - 1).fill (sp.getMass());

        for (size_t i = 0; i < Natom; i++)
        {
            create_electron(&sp, i);
        }

        DirectSum();

        double  tStatistics = 0.0;
        for (int i = 1; i <= ndt; i++, t = tstart + dt * i)
        {
            if (t>=tStatistics) {
                size_t  j;
                
                j = floor(tStatistics/histogram_dt);
                
                statTime(j) = tStatistics;
                
                for (size_t k=Natom; k<Nparticles; k++) {
                    double Ebody = 0.5*mass[k]*(vx[k]*vx[k] + vy[k]*vy[k] + vz[k]*vz[k]) + q[k]*phi[k] ;
                    if ( Ebody >=0 ) {
                        nfree_electrons[j] += 1.0;
                    }
                    if (revangle[k - Natom] > LocalizationAngle) {
                        nlocalized_electrons[j]++;
                    }
                }
                
                ntot_electrons(j) += Nparticles - Natom;
                
                cm(0, j) += mean(x.rows(0, Nparticles-1));
                cm(1, j) += mean(y.rows(0, Nparticles-1));
                cm(2, j) += mean(z.rows(0, Nparticles-1));
                
                vcm(0, j) += mean(vx.rows(0, Nparticles-1));
                vcm(1, j) += mean(vy.rows(0, Nparticles-1));
                vcm(2, j) += mean(vz.rows(0, Nparticles-1));
                
                Ekinavg(j) += kinetic_energy();
                Epotavg(j) += Utot;
                EoffsetAvg(j) += Eoffset;
                
                tStatistics += histogram_dt;
            }
            histogramNo = floor((i * dt) / histogram_dt);
            vec invmass = 1.0 / mass;

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
        
            Localization(dt);

            for (size_t j = 0; j < Natom; j++)
            {
                if (nloc[j] == 0)
                {
                    create_electron (&sp, j);
                }
            }

            if (i % ndtSnapshot == 0)
            {
                h5dump.newGroup();
                h5dump.add (fDumpVectors, Nparticles);
                h5dump.add ("nloc", nloc, Natom);
                h5dump.add ("next_atom", next_atom, Nparticles - Natom);
                h5dump.add ("revangle", revangle, Nparticles - Natom);
                h5dump.closeGroup();
            }
        }
    }
    
    quasi_free_hist *=  dt/(histogram_dt*(double)Natom*Nruns);
    valence_hist *=  dt/(histogram_dt*(double)Natom*Nruns);
    
    mat quasi_free_rho(quasi_free_hist), valence_rho(valence_hist);
    vec shellVolumes(histogram_bins);
    int Nbins = shellVolumes.n_rows;
    
    shellVolumes.rows(0, Nbins-2) = pow(histogram_bins.rows(1, Nbins-1), 3) -
        pow(histogram_bins.rows(0, Nbins-2), 3);
    shellVolumes(Nbins-1) = pow(histogram_rmax, 3) - pow(histogram_bins(Nbins-1), 3);
    shellVolumes *= 4.0*M_PI/3.0;
    
    for (size_t j=0; j<quasi_free_rho.n_cols; j++) {
        quasi_free_rho.col(j) = quasi_free_hist.col(j) / shellVolumes;
        valence_rho.col(j) = valence_hist.col(j) / shellVolumes;
    }
        

    
    
    h5dump.addStatistics("quasi_free_hist", quasi_free_hist);
    h5dump.addStatistics("quasi_free_rho", quasi_free_hist);
    h5dump.addStatistics("valence_hist", valence_hist);
    h5dump.addStatistics("valence_rho", valence_hist);
    h5dump.addStatistics("histogram_bins", histogram_bins);
    
    for (map<string, mat *>::iterator j=statFields.begin(); j != statFields.end(); j++) {
        *(j->second) /= (double)Nruns;
        h5dump.addStatistics(j->first.c_str(), *(j->second));
    }


    return 0;
}

