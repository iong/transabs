//
//  Atom.h
//  transabs
//
//  Created by Ionuț Georgescu on 8/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef ATOM_H
#define ATOM_H

class Atom
{
public:
    virtual double Eip ( int q ) = 0;
    virtual double getMass ( ) = 0;
};


#endif
