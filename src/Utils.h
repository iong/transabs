//
//  Polar.h
//  transabs
//
//  Created by Ionuț Georgescu on 8/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef transabs_Utils_h
#define transabs_Utils_h


#include <armadillo>
using namespace arma;

extern double  myrand(void);

extern void polar2cart ( double p[], double& x, double& y, double& z );

void    update_histogram(mat &h, size_t col, double rmax, double dr, vec &v);
//void    update_histogram(vec &h, double rmax, double dr, double v);

#endif
