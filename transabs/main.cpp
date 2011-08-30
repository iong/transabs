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
#include <boost/program_options.hpp>


#include "Eigen3.h"
#include "Argon.h"
#include "Utils.h"

using namespace std;
using namespace Eigen;
namespace po = boost::program_options;

po::variables_map vm;


int Natom;
Matrix3Xd r, v, f;
VectorXd    q, m;

double eps;
Argon   sp;


void process_options (int argc, char * argv[])
{
    try {

        po::options_description cmdline("Command line options");
        cmdline.add_options()
        ("help", "get help")
        ("config", po::value<string>(), "configuration file");

        po::options_description generic("General options");
        generic.add_options()
        ("pump.wave-length", po::value<double>()->required(), "wave length of the pump laser [nm]")
        ("pump.intensity", po::value<double>()->required(), "intensity of the pump laser [W/cm^2]")
        ("pump.duration", po::value<double>()->required(), "pulse length of the pump laser [fs]")
        ("coordinates", po::value<string>()->required(), "initial coordinates of the cluster atoms")
        ;

        cmdline.add(generic);

        po::positional_options_description  p;
        p.add("coordinates", 1);

        po::store(po::command_line_parser(argc, argv).options(cmdline).positional(p).run(), vm);

        if (vm.count("config")) {
            ifstream cfgin(vm["config"].as<string>().c_str());
            po::store(po::parse_config_file(cfgin, generic), vm);
            cfgin.close();
        }

        po::notify(vm);

        if (vm.count("help")) {
            cout << "\n\t\tAsk Ionut!\n\n"<<cmdline<<endl;
            exit(EXIT_FAILURE);
        }

        cout << vm["pump.duration"].as<double>()<<endl;
        }
        catch (po::required_option& ro) {
            cerr << "Required option " << ro.get_option_name() << " is missing.\n";
            exit(EXIT_FAILURE);
        }
        catch(exception& e) {
            cerr << "error: " << e.what() << "\n";
            exit(EXIT_FAILURE);
        }
}


int main (int argc, char * argv[])
{
    ifstream f;
    Matrix3Xd r(3, 4);

    process_options(argc, argv);

    f.open(argv[1]);
    f >> r;
    f.close();



    double pol[3];
    for (int i=0; i<Natom; i++) {
      Vector3d pol(sqrt(2.0 * (1.0/eps - sp.Eip(0))), acos(2.0*myrand() - 1.0),
                        2.0*M_PI*myrand() );
        v.col(i) = polar2cart(pol);

    }

    ndt = (tstart - tstop)/dt;
    for ( int i=0; i<=ndt; i++, t=tstart + dt*ndt) {
        v = v + (0.5*dt)* (m.asDiagonal().inverse()*f);
        r = r + dt*v;
        direct_force(q, r, U, f, Utot);
        v = v + (0.5*dt)* (m.asDiagonal().inverse()*f);
    }


    return 0;
}

