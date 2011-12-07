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

extern double t, tstart, tstop, dt;

extern vec    q, mass, invmass, x, y, z, vx, vy, vz, fx, fy, fz, phi, next_atom_dist,
    new_next_atom_dist, revangle;
extern ivec next_atom, new_next_atom, nlocByRevAngle, nlocByEnergy, valence;
extern double  Utot, Eoffset, rcluster;

extern void resize_vectors ();

extern int mpiNProcs, mpiRank;


#endif
