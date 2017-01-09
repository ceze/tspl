/*****************************************************************************
 *                             fastconv_usefftw_test.cpp
 *
 * Fast convolution testing by using FFTW.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <convolution_usefftw.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 3;
const   int     N = 5;


int main()
{
    Vector<Type> xn( M ), yn( N ), zn;

    for( int i=0; i<M; ++i )
        xn[i] = i;
    for( int i=0; i<N; ++i )
        yn[i] = i-N/2;

    // covolution
    zn = fastConvFFTW( xn, yn );
    cout << "xn:  " << xn << endl << endl
         << "yn:  " << yn << endl << endl;
    cout << "convolution:   " << zn << endl << endl;

    return 0;
}
