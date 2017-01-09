/*****************************************************************************
 *                                   rls_test.cpp
 *
 * RLS adaptive filter testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <convolution.h>
#include <vectormath.h>
#include <random.h>
#include <rls.h>


using namespace std;
using namespace splab;


typedef float   Type;
const   int     N = 1000;
const   int     orderRls = 1;
const   int     orderLrls = 16;
const   int     orderTrls = 12;
const   int     orderQrrls = 8;
const   int     sysLen = 8;
const   int     dispNumber = 10;


int main()
{
    int start = max(0,N-dispNumber);
    Vector<Type>    dn(N), xn(N), yn(N), sn(N), rn(N), en(N),
                    hn(sysLen+1), gn(orderLrls+1), wn(orderRls+1);
    Type lambda, delta, eps;


    cout << "/**************  Conventional RLS  <--->  Waveform Tracking \
*************/" << endl << endl;
    for( int k=0; k<N; ++k )
    {
        dn[k] = Type(sin(TWOPI*k/7));
        xn[k] = Type(cos(TWOPI*k/7));
    }
    lambda = Type(0.99);
    delta = dotProd(xn,xn)/N;

    cout << "The last " << dispNumber << " iterations result:" << endl << endl;
    cout << "observed" << "\t" << "desired" << "\t\t" << "output" << "\t\t"
         << "adaptive filter" << endl << endl;
    for( int k=0; k<start; ++k )
        yn[k] = rls( xn[k], dn[k], wn, lambda, delta );
    for( int k=start; k<N; ++k )
    {
        yn[k] = rls( xn[k], dn[k], wn, lambda, delta );
        cout << setiosflags(ios::fixed) << setprecision(4)
             << xn[k] << "\t\t" << dn[k] << "\t\t" << yn[k] << "\t\t";
        for( int i=0; i<=orderRls; ++i )
            cout << wn[i] << "\t";
        cout << endl;
    }
    cout << endl << "The theoretical optimal filter is:\t\t" << "-0.7972\t1.2788"
         << endl << endl << endl;


    cout << "/**************  Lattice RLS  <--->  Channel Equalization \
***************/" << endl << endl;
    for( int k=0; k<=sysLen; ++k )
        hn[k] = Type( 0.1 * pow(0.5,k) );
    dn = randn( 37, Type(0.0), Type(1.0), N );
    xn = wkeep( conv(dn,hn), N, "left" );
    lambda = Type(0.99), eps = Type(0.1);

    wn.resize(orderLrls);
    wn = Type(0.0);
    for( int k=0; k<N; ++k )
//        yn[k] = lrls( xn[k], dn[k], wn, lambda, eps, "on" );
        yn[k] = eflrls( xn[k], dn[k], wn, lambda, eps, "on" );
    Vector<Type> Delta(orderLrls+1);
    Delta[(orderLrls+1)/2] = Type(1.0);
    for( int k=0; k<=orderLrls; ++k )
//        gn[k] = lrls( Delta[k], Type(0.0), wn, lambda, eps, "off" );
        gn[k] = eflrls( Delta[k], Type(0.0), wn, lambda, eps, "off" );
    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "The original system:   " << hn << endl;
    cout << "The inverse system:   " << gn << endl;
    cout << "The cascade system:   " << conv( gn, hn ) << endl << endl;
//

    cout << "/************  Transversal RLS  <--->  System Identification \
************/" << endl << endl;
    Vector<Type> sys(8);
    sys[0] = Type(0.1);   sys[1] = Type(0.3);
    sys[2] = Type(0.0);   sys[3] = Type(-0.2);
    sys[4] = Type(-0.4);  sys[5] = Type(-0.7);
    sys[6] = Type(-0.4);  sys[7] = Type(-0.2);
    xn = randn( 37, Type(0.0), Type(1.0), N );
    dn = wkeep( conv(xn,sys), N, "left" );
    lambda = Type(0.99), eps = Type(1.0);
    wn.resize(orderTrls);
    wn = Type(0.0);

    for( int k=0; k<start; ++k )
        yn[k] = sftrls( xn[k], dn[k], wn, lambda, eps, "on" );
    cout << "The last " << dispNumber << " iterations result:" << endl << endl;
    cout << "input signal" << "    " << "original system output" << "    "
         << "identified system output" << endl << endl;
    for( int k=start; k<N; ++k )
    {
        yn[k] = sftrls( xn[k], Type(0.0), wn, lambda, eps, "off" );
        cout << setiosflags(ios::fixed) << setprecision(4)
             << xn[k] << "\t\t\t" << dn[k] << "\t\t\t" << yn[k] << endl;
    }
    cout << endl << "The unit impulse response of original system:   "
         << sys << endl;
    cout << "The unit impulse response of identified system:   "
         << wn << endl << endl;


    cout << "/***************  QR Based RLS  <--->  Signal Enhancement \
 ***************/" << endl << endl;
    sn = sin( linspace( Type(0.0), Type(4*TWOPI), N ) );
    rn = randn( 37, Type(0.0), Type(1.0), N );
    dn = sn + rn;
    int delay = orderQrrls/2;
    for( int i=0; i<N-delay; ++i )
        xn[i] = rn[i+delay];
    for( int i=N-delay; i<N; ++i )
        xn[i] = 0;
    Type lambdaSqrt = Type(1.0);
    wn.resize(orderQrrls);
    wn = Type(0.0);

    for( int k=0; k<N; ++k )
        yn[k] = qrrls( xn[k], dn[k], wn, lambdaSqrt, "on" );
    en = dn - yn;
    for( int i=0; i<N/2; ++i )
    {
        sn[i] = 0;
        rn[i] = 0;
        en[i] = 0;
    }
    cout << "The last " << dispNumber << " iterations result:" << endl << endl;
    cout << "noised signal\t\t" << "enhanced signal\t\t" << "original signal"
         << endl << endl;
    for( int k=start; k<N; ++k )
        cout << setiosflags(ios::fixed) << setprecision(4)
             << dn[k] << "\t\t\t" << en[k] << "\t\t\t" << sn[k] << endl;
    cout << endl << "The SNR before denoising is:   "
         << 20*(log10(norm(sn))/log10(norm(rn))) << " dB" << endl;
    cout << "The SNR after denoising is:    "
         << 20*(log10(norm(en))/log10(norm(sn-en))) << " dB" << endl << endl;


    return 0;
}
