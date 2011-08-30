//
//  DirectSumation.cpp
//  transabs
//
//  Created by Ionuț Georgescu on 8/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "Eigen3.h"

using namespace Eigen;

void DirectSumation(VectorXd q, Matrix3Xd& r, VectorXd& pot, Matrix3Xd& f,
                    double& Utot)
{
  int i, j, N;
  
  N = r.cols();
  for (i=0; i<N; i++) {
    dr.leftCols(i).colwise() = r.leftCols(i).colwise - r.col(i);
    rsq.head(i) = dr.leftCols(i).square().colwise().sum();
    rabs.head(i) = (rsq.head(i) + eps2).sqrt();
    df.leftCols(i) = (q/(rsq*rqbs)).asDiagonal() * rsqdr.leftCols(i);
  }
}