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


void    update_histogram(mat &h, size_t col, double rmax, double dr, vec &v)
{
    size_t  i;
    
    for (size_t j=0; j<v.n_rows; j++) {
        if (v[j]>=rmax) continue;
        i = floor(v[j]/dr);
        /*
        if (v[j]>=rmax) {
            i = h.n_rows-1;
        }
         */
        h(i, col) += 1.0;
    }
}

void polar2cart (double p[], double& x, double& y, double& z)
{
        vec dest;

        x = ( p[0] * sin ( p[1] ) ) * cos ( p[2] );
        y = ( p[0] * sin ( p[1] ) ) * sin ( p[2] );
        z = p[0] * cos ( p[1] );
}