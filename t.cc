#include <iostream>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

int main ( int argc, char *argv[] )
{
    int i, j;
    Array3Xd f ( 3, 5 );
    Array3Xd v ( 3, 5 );
    ArrayXd m ( 5 );

    for ( i = 0; i < 5; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            f ( j, i ) = i * j;
        }
    }

    m << 1, 2, 1, 2, 1;


    for ( i = 0; i < f.cols(); i++ ) {
        v.col ( i ) = f.col ( i ) / m ( i );
    }
    cout << v << endl;

    f.resize ( NoChange, 2 );
    //m.resize(3);
}
