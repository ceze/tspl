/*****************************************************************************
 *                              convolution.cpp
 *
 * Convolution testing.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <convolution.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 3;
const   int     N = 5;


int main()
{
    Vector<Type> xn( M ), yn( N );
    Vector<Type> zn;

    for( int i=0; i<M; ++i )
        xn[i] = i;
    for( int i=0; i<N; ++i )
        yn[i] = i-N/2;

    // convolution
    zn = conv( xn, yn );
    cout << "xn:  " << xn << endl << "yn:  " << yn << endl;
    cout << "convolution of xn and yn:   " << zn << endl;
    zn = fastConv( xn, yn );
    cout << "fast convolution of xn and yn:   " << zn << endl;

    return 0;
}
