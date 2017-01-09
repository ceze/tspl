/*****************************************************************************
 *                                  lms_test.cpp
 *
 * LMS adaptive filter testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <lms.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 50;
const   int     order = 1;
const   int     dispNumber = 10;


int main()
{
    Vector<Type> dn(N), xn(N), yn(N), wn(order+1);
    for( int k=0; k<N; ++k )
    {
        xn[k] = Type(cos(TWOPI*k/7));
        dn[k] = Type(sin(TWOPI*k/7));
    }
    int start = max(0,N-dispNumber);
    Type xnPow = dotProd(xn,xn)/N;
    Type mu = Type( 0.1 / ((order+1)*xnPow) );
    Type rho = Type(1.0), gamma = Type(1.0e-9), alpha = Type(0.08);

    cout << "The last " << dispNumber
         << " iterations of Conventional-LMS:" << endl;
    cout << "observed" << "\t" << "desired" << "\t\t" << "output" << "\t\t"
         << "adaptive filter" << endl << endl;
    wn = 0;
    for( int k=0; k<start; ++k )
        yn[k] = lms( xn[k], dn[k], wn, mu );
    for( int k=start; k<N; ++k )
    {
        yn[k] = lms( xn[k], dn[k], wn, mu );
        cout << setiosflags(ios::fixed) << setprecision(4)
             << xn[k] << "\t\t" << dn[k] << "\t\t" << yn[k] << "\t\t";
        for( int i=0; i<=order; ++i )
            cout << wn[i] << "\t";
        cout << endl;
    }
    cout << endl << endl;

    cout << "The last " << dispNumber << " iterations of LMS-Newton:" << endl;
    cout << "observed" << "\t" << "desired" << "\t\t" << "output" << "\t\t"
         << "adaptive filter" << endl << endl;
    wn = 0;
    for( int k=0; k<start; ++k )
        yn[k] = lmsNewton( xn[k], dn[k], wn, mu, alpha, xnPow );
    for( int k=start; k<N; ++k )
    {
        yn[k] = lmsNewton( xn[k], dn[k], wn, mu, alpha, xnPow );
        cout << setiosflags(ios::fixed) << setprecision(4)
             << xn[k] << "\t\t" << dn[k] << "\t\t" << yn[k] << "\t\t";
        for( int i=0; i<=order; ++i )
            cout << wn[i] << "\t";
        cout << endl;
    }
    cout << endl << endl;

    cout << "The last " << dispNumber
         << " iterations of Normalized-LMS:" << endl;
    cout << "observed" << "\t" << "desired" << "\t\t" << "output" << "\t\t"
         << "adaptive filter" << endl << endl;
    wn = 0;
    for( int k=0; k<start; ++k )
        yn[k] = lmsNormalize( xn[k], dn[k], wn, rho, gamma );
    for( int k=start; k<N; ++k )
    {
        yn[k] = lmsNormalize( xn[k], dn[k], wn, rho, gamma );
        cout << setiosflags(ios::fixed) << setprecision(4)
             << xn[k] << "\t\t" << dn[k] << "\t\t" << yn[k] << "\t\t";
        for( int i=0; i<=order; ++i )
            cout << wn[i] << "\t";
        cout << endl;
    }
    cout << endl << endl;

    cout << "The theoretical optimal filter is:\t\t" << "-0.7972\t1.2788"
         << endl << endl;

    return 0;
}
