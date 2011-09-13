#ifndef transabs_Common_h
#define transabs_Common_h

#include <string>

#include <boost/program_options.hpp>

#include <armadillo>



using namespace std;
using namespace arma;
namespace po = boost::program_options;

extern po::variables_map vm;


extern double soft_core, LocalizationAngle;

extern size_t Natom, Nparticles;
extern vec    q, mass, invmass, x, y, z, vx, vy, vz, fx, fy, fz, phi, next_atom_dist,
    new_next_atom_dist, revangle, next_atom_t0;
extern ivec next_atom, new_next_atom, nloc, valence;
extern double  Utot, Eoffset;

extern mat quasi_free_hist, valence_hist;
extern double histogram_rmax;
extern double histogram_dr;
extern size_t histogramNo;
extern vec histogram_norm;

extern void resize_vectors ();


#endif