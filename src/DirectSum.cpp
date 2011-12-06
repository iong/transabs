//
//  DirectSumation.cpp
//  transabs
//
//  Created by Ionuț Georgescu on 8/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "Common.h"
#include "DirectSum.h"
#include "Utils.h"
#include "Statistics.h"

void DirectSum ()
{
    size_t i;
    u32 closest_atom_id;
    vec    dphi (1);
    double  eps2, closest_atom, risq, rcutoffsq = rcluster * rcluster * 4.0;
    vec dx (Nparticles);
    vec dy (Nparticles);
    vec dz (Nparticles);
    vec ir (Nparticles);
    vec ir2 (Nparticles);
    vec ir3 (Nparticles);
    vec atomDist (Natom);

    fx.zeros();
    fy.zeros();
    fz.zeros();
    phi.zeros();
    Utot = 0.0;

    eps2 = soft_core * soft_core;

    new_next_atom_dist.fill (INFINITY);
    vec realCharge = q.subvec(allAtoms) - 1.0;
    vec atomIDs=linspace<vec>(0, Natom-1, Natom);
    
    for (i = Nparticles - 1; i > 0; i--)
    {
        dx.rows(0, i-1) = x.rows(0, i-1) - x (i);
        dy.rows(0, i-1) = y.rows(0, i-1) - y (i);
        dz.rows(0, i-1) = z.rows(0, i-1) - z (i);

        ir2.rows(0, i-1) = square (dx.rows(0, i-1)) + square (dy.rows(0, i-1)) + square (dz.rows(0, i-1));

        risq = x (i) * x (i) + y (i) * y (i) + z (i) * z (i);

        if (i >= Natom && risq < rcutoffsq)
        {
            atomDist = sqrt (ir2.rows (0, Natom - 1));

            incrementRadialDistributionsQ(i, atomDist, realCharge);

            closest_atom = atomDist.min (closest_atom_id);

            if (closest_atom < new_next_atom_dist[i])
            {
                new_next_atom_dist[i] = closest_atom;
                new_next_atom[i] = closest_atom_id;
            }
        }
        else if (i >= Natom)
        {
            new_next_atom[i] = -1;
        }

        ir2.rows(0, i-1) = 1.0 / (ir2.rows(0, i-1) + eps2);

        ir.rows(0, i-1) = sqrt (ir2.rows(0, i-1));
        ir3.rows(0, i-1) = ir.rows(0, i-1) % ir2.rows(0, i-1);

        phi.rows(0, i-1) += q (i) * ir.rows(0, i-1);

        dphi = sum (ir.rows(0, i-1) % q.rows(0, i-1));

        phi (i) += dphi (0);
        Utot += dphi (0) * q (i);

        dx.rows(0, i-1) = dx.rows(0, i-1) % ir3.rows(0, i-1);
        dy.rows(0, i-1) = dy.rows(0, i-1) % ir3.rows(0, i-1);
        dz.rows(0, i-1) = dz.rows(0, i-1) % ir3.rows(0, i-1);

        fx.rows(0, i-1) += dx.rows(0, i-1) * q (i) ;
        fy.rows(0, i-1) += dy.rows(0, i-1) * q (i) ;
        fz.rows(0, i-1) += dz.rows(0, i-1) * q (i) ;

        fx.rows (i, i) += - sum (dx.rows(0, i-1) % q.rows(0, i-1));
        fy.rows (i, i) += - sum (dy.rows(0, i-1) % q.rows(0, i-1));
        fz.rows (i, i) += - sum (dz.rows(0, i-1) % q.rows(0, i-1));
    }
}
