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
#include "DirectSum.h"
#include "Utils.h"

using namespace std;

po::variables_map vm;

int Natom;
vec    q, mass, invmass, x, y, z, vx, vy, vz, fx, fy, fz, phi;
double  Utot;

Argon   sp;

double t, tstart, tstop, dt;
double soft_core;


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
            ( "pump.duration", po::value<double>()->required(), "pulse length of the pump laser [fs]" )
            ( "coordinates", po::value<string>()->required(), "initial coordinates of the cluster atoms" )
            ( "eps", po::value<double>()->required(), "soft-core parameter")
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
}


int main ( int argc, char * argv[] )
{
    mat rt;
    int ndt;

    process_options ( argc, argv );

    rt.load ( vm["coordinates"].as<string>() );
    Natom = rt.n_rows;

    x = rt.col(0); y = rt.col(1); z = rt.col(2);
    vx.zeros (Natom); vy = vx; vz = vx;
    fx=vx; fy = vx; fz = vx;
    q.ones ( Natom );
    mass.ones ( Natom );

    invmass = 1.0 / mass;

    rt.reset();

    double pol[3];
    for ( int i = 0; i < Natom; i++ ) {
        pol [0] = sqrt ( 2.0 * ( 1.0 / soft_core - sp.Eip ( 0 ) ) );
        pol [1] = acos ( 2.0 * myrand() - 1.0 );
        pol [2] = 2.0 * M_PI * myrand();
        polar2cart ( pol, vx(i), vy(i), vz(i) );

    }

    ndt = ( tstart - tstop ) / dt;
    for ( int i = 0; i <= ndt; i++, t = tstart + dt * ndt ) {
        vx += ( 0.5 * dt ) * fx % invmass;
        vy += ( 0.5 * dt ) * fy % invmass;
        vz += ( 0.5 * dt ) * fz % invmass;

        x += dt * vx;
        y += dt * vy;
        z += dt * vz;

        DirectSum(q, x, y, z, fx, fy, fz, phi, Utot);

        vx += ( 0.5 * dt ) * fx % invmass;
        vy += ( 0.5 * dt ) * fy % invmass;
        vz += ( 0.5 * dt ) * fz % invmass;
    }


    return 0;
}

