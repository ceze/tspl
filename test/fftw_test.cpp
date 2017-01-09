/*****************************************************************************
 *                               fftw_test.cpp
 *
 * FFTW interface testing.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <fftw.h>


using namespace std;
using namespace splab;


typedef float   Type;
const   int     N = 7;


int main()
{
    // complex to complex dft test...
    Vector< complex<Type> > xn( N );
    Vector< complex<Type> > yn( N );
    Vector< complex<Type> > Xk( N );

    for( int i=0; i<N; ++i )
    {
        Type theta = Type( 2*PI * i / N );
        xn[i] = complex<Type>( cos(theta), sin(theta) );
    }

    cout << setiosflags(ios::fixed) << setiosflags(ios::showpos) << setprecision(8);
    cout << "xn:   " << xn << endl << endl;

    fftw( xn, Xk );
    cout << "Xk=fft(xn):   " << Xk << endl << endl;
    ifftw( Xk, yn );
    cout << "xn-ifft(Xk):   " << xn-yn << endl << endl;

    // real to complex and complex to real dft test...
    Vector<Type> sn(N), tn(N);
    Vector< complex<Type> > Sk;
    for( int i=0; i<N; ++i )
    {
        Type theta = Type( 2*PI * i / N );
        sn[i] = sin(theta);
    }
    cout << "sn:   " << sn << endl;

    fftw( sn, Sk );
    cout << "Sk=fft(sn):   " << Sk << endl << endl;
    ifftw( Sk, tn );
    cout << "sn-ifft(Sk):   " << sn-tn << endl;

    return 0;
}
