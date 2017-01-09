/*****************************************************************************
 *                                   kalman.h
 *
 * Kalman filter testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <kalman.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 2;
const   int     M = 2;
const   int     T = 20;


int main()
{
    Matrix<Type> A(N,N), C(M,N), Q(N,N), R(M,M);
    A = eye( N, Type(1.0) );    C = eye( N, Type(1.0) );
    Q = eye( N, Type(1.0) );    R = eye( N, Type(2.0) );

    Vector<Type> x(N,Type(1.0)), y(M), ytInit(M);
    ytInit[0] = Type(0.5);  ytInit[1] = Type(2.0);
    Matrix<Type> yt(M,T);
    for( int t=0; t<T; ++t )
        yt.setColumn( ytInit, t );

    Vector<Type> intV( N, Type(10.0) );
    for( int t=0; t<T; ++t )
    {
        y = yt.getColumn(t);
        x = kalman( A, C, Q, R, x, y, intV );
        cout << "Estimation of xt at the " << t << "th iteratin:   " << x << endl;
    }

    cout << "The theoretical xt should converge to:   " << ytInit << endl;

    return 0;
}
