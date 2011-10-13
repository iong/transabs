//
//  Atom.h
//  transabs
//
//  Created by Ionuț Georgescu on 8/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef ATOM_H
#define ATOM_H

/**
 * J. Chem. Phys. 118, 4976 (2003); doi:10.1063/1.1543944
 */
class Atom
{
public:
    virtual double Eip ( int q ) = 0;
    virtual double getMass ( ) = 0;
    virtual double getLJDimer() = 0;
};


#endif
