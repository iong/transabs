//
//  DirectSumation.cpp
//  transabs
//
//  Created by Ionuț Georgescu on 8/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//



#include "Common.h"
#include "DirectSum.h"

void DirectSum ( vec& q, vec& x, vec& y, vec& z, vec& fx, vec& fy, vec& fz,
                 vec& phi, double& Utot )
{
    int i, j, N;
    vec    dphi(1);
    double  eps2;

    N = q.n_rows;

    vec dx ( N ), dy ( N ), dz ( N ), ir ( N ), ir2 ( N ), ir3 ( N );



    fx = 0.0;
    fy = 0.0;
    fz = 0.0;
    phi = 0.0;
    Utot = 0.0;

    eps2 = soft_core * soft_core;

    for ( i = N - 1; i > 0; i-- ) {
        dx.rows ( 0, i - 1 ) = x.rows ( 0, i - 1 ) - x ( i );
        dy.rows ( 0, i - 1 ) = y.rows ( 0, i - 1 ) - y ( i );
        dz.rows ( 0, i - 1 ) = z.rows ( 0, i - 1 ) - z ( i );

        ir2.rows ( 0, i - 1 ) = 1.0 / (
                    ( dx.rows ( 0, i - 1 ) % dx.rows ( 0, i - 1 ) )
                    + ( dy.rows ( 0, i - 1 ) % dy.rows ( 0, i - 1 ) )
                    + ( dz.rows ( 0, i - 1 ) % dz.rows ( 0, i - 1 ) )
                    + eps2 );

        ir.rows ( 0, i - 1 ) = sqrt ( ir2.rows ( 0, i - 1 ) );
        ir3.rows ( 0, i - 1 ) = ir.rows ( 0, i - 1 ) % ir2.rows ( 0, i - 1 );

        phi.rows ( 0, i - 1 ) += q ( i ) * ir.rows ( 0, i - 1 );

        dphi = sum ( ir.rows ( 0, i - 1 ) % q.rows ( 0, i - 1 ) );

        phi ( i ) += dphi(0);
        Utot += dphi(0) * q ( i );

        dx.rows ( 0, i - 1 ) = dx.rows ( 0, i - 1 ) % ir3.rows ( 0, i - 1 );
        dy.rows ( 0, i - 1 ) = dy.rows ( 0, i - 1 ) % ir3.rows ( 0, i - 1 );
        dz.rows ( 0, i - 1 ) = dz.rows ( 0, i - 1 ) % ir3.rows ( 0, i - 1 );

        fx.rows ( 0, i - 1 ) += dx.rows ( 0, i - 1 ) * q ( i ) ;
        fy.rows ( 0, i - 1 ) += dy.rows ( 0, i - 1 ) * q ( i ) ;
        fz.rows ( 0, i - 1 ) += dz.rows ( 0, i - 1 ) * q ( i ) ;

        fx.rows ( i,i ) = - sum ( dx.rows ( 0, i - 1 ) % q.rows ( 0, i - 1 ) );
        fy.rows ( i,i ) = - sum ( dy.rows ( 0, i - 1 ) % q.rows ( 0, i - 1 ) );
        fz.rows ( i,i ) = - sum ( dz.rows ( 0, i - 1 ) % q.rows ( 0, i - 1 ) );
    }
}
