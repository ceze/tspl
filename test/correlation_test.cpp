/*****************************************************************************
 *                             correlation_test.cpp
 *
 * Correlation testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <correlation.h>


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

    // auto and cross correlation functions
    cout << "auto-correlation of xn:   " << corr(xn) << endl;
    cout << "biased auto-correlation of xn:   " << corr(xn,"biased") << endl;
    cout << "unbiased auto-correlation of xn:   " << corr(xn,"unbiased") << endl;

    cout << "cross-correlation of xn and yn:   " << corr(xn,yn) << endl;
    cout << "cross-correlation of yn and xn:   " << corr(yn,xn) << endl;

    cout << "biased cross-correlation of xn and yn:   "
         << corr(xn,yn,"biased") << endl;
    cout << "biased cross-correlation of yn and xn:   "
         << corr(yn,xn,"biased") << endl;

    cout << "unbiased cross-correlation of xn and yn:   "
         << corr(xn,yn,"unbiased") << endl;
    cout << "unbiased cross-correlation of yn and xn:   "
         << corr(yn,xn,"unbiased") << endl;

    // fast auto and cross correlation functions
    cout << "fast auto-correlation of xn:   "
         << fastCorr(xn) << endl;
    cout << "fast biased auto-correlation of xn:   "
         << fastCorr(xn,"biased") << endl;
    cout << "fast unbiased auto-correlation of xn:   "
         << fastCorr(xn,"unbiased") << endl;

    cout << "fast cross-correlation of xn and yn:   " << fastCorr(xn,yn) << endl;
    cout << "fast cross-correlation of yn and xn:   " << fastCorr(yn,xn) << endl;

    cout << "fast biased cross-correlation of xn and yn:   "
         << fastCorr(xn,yn,"biased") << endl;
    cout << "fast biased cross-correlation of yn and xn:   "
         << fastCorr(yn,xn,"biased") << endl;

    cout << "fast unbiased cross-correlation of xn and yn:   "
         << fastCorr(xn,yn,"unbiased") << endl;
    cout << "fast unbiased cross-correlation of yn and xn:   "
         << fastCorr(yn,xn,"unbiased") << endl;

    return 0;
}
