//
//  Utils.cpp
//  transabs
//
//  Created by Ionuț Georgescu on 8/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


#include <cmath>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include "Utils.h"

static boost::mt19937 rng;

double  myrand(void)
{
    boost::uniform_01<double> rnd01;

    return rnd01(rng);
}


void polar2cart (double p[], double& x, double& y, double& z)
{
        vec dest;

        x = ( p[0] * sin ( p[1] ) ) * cos ( p[2] );
        y = ( p[0] * sin ( p[1] ) ) * sin ( p[2] );
        z = p[0] * cos ( p[1] );
}