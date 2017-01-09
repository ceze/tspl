/*****************************************************************************
 *                             fastcorr_usefftw_test.cpp
 *
 * Fast correlation by using FFTW testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <correlation_usefftw.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 3;
const   int     N = 5;


int main()
{
    Vector<Type> xn( M ), yn( N );

    for( int i=0; i<M; ++i )
        xn[i] = i;
    for( int i=0; i<N; ++i )
        yn[i] = i-N/2;

    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "xn:   " << xn << endl;
    cout << "yn:   " << yn << endl;

    // fast auto and cross correlation functions by using FFTW
    cout << "fast auto-correlation of xn:   "
         << fastCorrFFTW(xn) << endl;
    cout << "fast biased auto-correlation of xn:   "
         << fastCorrFFTW(xn,"biased") << endl;
    cout << "fast unbiased auto-correlation of xn:   "
         << fastCorrFFTW(xn,"unbiased") << endl;

    cout << "fast cross-correlation of xn and yn:   "
         << fastCorrFFTW(xn,yn) << endl;
    cout << "fast cross-correlation of yn and xn:   "
         << fastCorrFFTW(yn,xn) << endl;

    cout << "fast biased cross-correlation of xn and yn:   "
         << fastCorrFFTW(xn,yn,"biased") << endl;
    cout << "fast biased cross-correlation of yn and xn:   "
         << fastCorrFFTW(yn,xn,"biased") << endl;

    cout << "fast unbiased cross-correlation of xn and yn:   "
         << fastCorrFFTW(xn,yn,"unbiased") << endl;
    cout << "fast unbiased cross-correlation of yn and xn:   "
         << fastCorrFFTW(yn,xn,"unbiased") << endl;

    return 0;
}
