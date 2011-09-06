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

#include "Common.h"
#include "Argon.h"
#include "Laser.h"
#include "DirectSum.h"
#include "Utils.h"

using namespace std;

po::variables_map vm;

int Natom, Nparticles, NMaxParticles;
vec    q, mass, invmass, x, y, z, vx, vy, vz, fx, fy, fz, phi;
double  Utot;

Argon   sp;

double t, tstart, tstop, dt;
double soft_core;

GaussianLaserPulse pump;


void process_options ( int argc, char * argv[] )
{
    try {

        po::options_description cmdline ( "Command line options" );
        cmdline.add_options()
        ( "help", "get help" )
        ( "config", po::value<string>(), "configuration file" );

        po::options_description generic ( "General options" );
        generic.add_options()
        ( "pump.wave-length", po::value<double>()->required(), "wave length of the pump laser [nm]" )
        ( "pump.intensity", po::value<double>()->required(), "intensity of the pump laser [W/cm^2]" )
        ( "pump.fwhm", po::value<double>()->required(), "pulse length of the pump laser [fs]" )
        ( "coordinates", po::value<string>()->required(), "initial coordinates of the cluster atoms" )
        ( "eps", po::value<double>()->required(), "soft-core parameter" )
        ( "tstart", po::value<double>()->required(), "soft-core parameter" )
        ( "tstop", po::value<double>()->required(), "soft-core parameter" )
        ( "dt", po::value<double>()->required(), "soft-core parameter" )
        ;

        cmdline.add ( generic );

        po::positional_options_description  p;
        p.add ( "coordinates", 1 );

        po::store ( po::command_line_parser ( argc, argv ).options ( cmdline ).positional ( p ).run(), vm );

        if ( vm.count ( "config" ) ) {
            ifstream cfgin ( vm["config"].as<string>().c_str() );
            po::store ( po::parse_config_file ( cfgin, generic ), vm );
            cfgin.close();
        }

        po::notify ( vm );

        if ( vm.count ( "help" ) ) {
            cout << "\n\t\tAsk Ionut!\n\n" << cmdline << endl;
            exit ( EXIT_FAILURE );
        }
    }
    catch ( po::required_option& ro ) {
        cerr << "Required option " << ro.get_option_name() << " is missing.\n";
        exit ( EXIT_FAILURE );
    }
    catch ( exception& e ) {
        cerr << "error: " << e.what() << "\n";
        exit ( EXIT_FAILURE );
    }
    soft_core = vm["eps"].as<double>();
    tstart = vm["tstart"].as<double>();
    tstop = vm["tstop"].as<double>();
    dt = vm["dt"].as<double>();

    // I_au = c eps_0 E_au^2/2
    // w_au = 2*pi*\hbar*c/E_H
    pump = GaussianLaserPulse ( vm["pump.intensity"].as<double>() / 3.50944493695457e+16,
                45.5633525101396 / vm["pump.wave-length"].as<double>(),
                vm["pump.fwhm"].as<double>() * 41.3413733524035 );
    pump.setZero ( 2.5 * vm["pump.fwhm"].as<double>() * 41.3413733524035 );
}

void DumpSnapshot(int   i)
{
    stringstream fn;
}

double kinetic_energy()
{
    vec ret ( 1 );

    ret.rows ( 0, 0 ) =  0.5 * sum (
                             ( square(vx.rows ( 0, Nparticles - 1 ) )
                               + square(vy.rows ( 0, Nparticles - 1 ))
                               + square(vz.rows ( 0, Nparticles - 1 )) ) % mass.rows ( 0, Nparticles - 1 ) ) ;

    return ret ( 0 );
}


int main ( int argc, char * argv[] )
{
    mat rt;
    int ndt;

    process_options ( argc, argv );

    rt.load ( vm["coordinates"].as<string>() );
    rt *= 3.84 * 1.89;
    Natom = rt.n_rows;
    NMaxParticles = 5 * Natom;

    x.zeros ( NMaxParticles );
    y = x;
    z = x;
    vx = x;
    vy = x;
    vz = x;
    fx = x;
    fy = x;
    fz = x;
    phi = x;
    q = x;
    mass = x;
    invmass = x;

    x.rows ( 0, Natom - 1 ) = rt.col ( 0 );
    y.rows ( 0, Natom - 1 ) = rt.col ( 1 );
    z.rows ( 0, Natom - 1 ) = rt.col ( 2 );

    q.rows ( 0, Natom - 1 ).fill ( 1.0 );
    mass.rows ( 0, Natom - 1 ).fill ( sp.getMass() );


    rt.reset();

    q.rows ( Natom, 2*Natom - 1 ).fill ( -1.0 );
    mass.rows ( Natom, 2*Natom - 1 ).fill ( 1.0 );
    x.rows ( Natom, 2*Natom - 1 ) = x.rows ( 0, Natom - 1 );
    y.rows ( Natom, 2*Natom - 1 ) = y.rows ( 0, Natom - 1 );
    z.rows ( Natom, 2*Natom - 1 ) = z.rows ( 0, Natom - 1 );

    Nparticles = 2 * Natom;

    double pol[3];
    for ( int i = 0; i < Natom; i++ ) {
        pol [0] = sqrt ( 2.0 * ( 1.0 / soft_core - sp.Eip ( 0 ) ) );
        pol [1] = acos ( 2.0 * myrand() - 1.0 );
        pol [2] = 2.0 * M_PI * myrand();
        polar2cart ( pol, vx ( Natom + i ), vy ( Natom + i ), vz ( Natom + i ) );

    }

    invmass.rows ( 0, Nparticles - 1 ) = 1.0 / mass.rows ( 0, Nparticles - 1 );

    double Ekin;
    ofstream eout ( "energy.dat" );

    DirectSum ();
    Ekin = kinetic_energy();
    eout << t << "\t" << Utot << "\t" << Ekin << "\t" << Utot + Ekin << endl;

    ndt = ( tstop - tstart ) / dt;
    for ( int i = 1; i <= ndt; i++, t = tstart + dt * i ) {
        vx.rows ( 0, Nparticles - 1 ) += ( 0.5 * dt ) * fx.rows ( 0, Nparticles - 1 ) % q.rows ( 0, Nparticles - 1 ) % invmass.rows ( 0, Nparticles - 1 );
        vy.rows ( 0, Nparticles - 1 ) += ( 0.5 * dt ) * fy.rows ( 0, Nparticles - 1 ) % q.rows ( 0, Nparticles - 1 ) % invmass.rows ( 0, Nparticles - 1 );
        vz.rows ( 0, Nparticles - 1 ) += ( 0.5 * dt ) * ( fz.rows ( 0, Nparticles - 1 ) + pump.field ( t ) ) % q.rows ( 0, Nparticles - 1 ) % invmass.rows ( 0, Nparticles - 1 );

        x.rows ( 0, Nparticles - 1 ) += dt * vx.rows ( 0, Nparticles - 1 );
        y.rows ( 0, Nparticles - 1 ) += dt * vy.rows ( 0, Nparticles - 1 );
        z.rows ( 0, Nparticles - 1 ) += dt * vz.rows ( 0, Nparticles - 1 );

        DirectSum ();


        vx.rows ( 0, Nparticles - 1 ) += ( 0.5 * dt ) * fx.rows ( 0, Nparticles - 1 ) % q.rows ( 0, Nparticles - 1 ) % invmass.rows ( 0, Nparticles - 1 );
        vy.rows ( 0, Nparticles - 1 ) += ( 0.5 * dt ) * fy.rows ( 0, Nparticles - 1 ) % q.rows ( 0, Nparticles - 1 ) % invmass.rows ( 0, Nparticles - 1 );
        vz.rows ( 0, Nparticles - 1 ) += ( 0.5 * dt ) * ( fz.rows ( 0, Nparticles - 1 ) + pump.field ( t ) ) % q.rows ( 0, Nparticles - 1 ) % invmass.rows ( 0, Nparticles - 1 );

        Ekin = kinetic_energy();
        eout << t << "\t"<< pump.field(t)<<"\t"<< Utot << "\t" << Ekin << "\t" << Utot + Ekin << endl;
    }
    eout.close();

    return 0;
}

