//
//  DenseBaseIStreamInitializer.h
//  transabs
//
//  Created by Ionuț Georgescu on 8/26/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef transabs_DenseBaseExtractOperator_h
#define transabs_DenseBaseExtractOperator_h


friend std::istream& operator>>(std::istream &sin, DenseBase<Derived> &m)
{
    Index r, c;
    
    r=0; c=0;
    for (r=0; r<m.rows(); r++) {
        for (c=0; c<m.cols(); c++) {
            sin >> m(r,c);
        }
    }
    return sin;
};


#endif
