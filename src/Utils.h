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

extern void init_rng(int skip);
extern double  myrand(void);

extern void polar2cart ( double p[], double& x, double& y, double& z );

#endif
