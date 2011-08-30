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

#include "Eigen3.h"

boost::mt19937 rng;

double  myrand(void)
{
    boost::uniform_01<double> rnd01;

    return rnd01(rng);
}


Eigen::Vector3d polar2cart ( Eigen::Vector3d& src )
{
        Eigen::Vector3d dest;

        dest(0) = ( src(0) * sin ( src(1) ) ) * cos ( src(2) );
        dest(1) = ( src(0) * sin ( src(1) ) ) * sin ( src(2) );
        dest(2) = src(0) * cos ( src(1) );
}