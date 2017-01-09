/*****************************************************************************
 *                                wiener_test.h
 *
 * Wiener filter testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <wiener.h>
#include <random.h>
#include <vectormath.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 1024;
const   int     M = 12;
const   int     fOrder = 8;
const   int     pOrder = 3;


int main()
{
    Vector<Type> tn(N), dn(N), vn(N), xn(N), yn(N);
    tn = linspace( Type(0.0), Type(2*TWOPI), N );
    dn = sin(tn);
    vn = randn( 37, Type(0.0), Type(1.0), N );
    xn = dn + vn;

    Vector<Type> hn = wienerFilter( xn, dn, fOrder );
    cout << "Wiener filter hn:   " << hn << endl;
    yn = wkeep( conv(xn,hn), N, "left" );
    cout << "original relative error:   " << norm(dn-xn)/norm(dn) << endl;
    cout << "filtered relative error:   " << norm(dn-yn)/norm(dn) << endl << endl;

    Vector<Type> sn(M);
    for( int i=0; i<M; ++i )
        sn[i] = sin( Type(i*TWOPI/10) );
    Vector<Type> pn = wienerPredictor( sn, pOrder );
    Vector<Type> snPred = wkeep( conv(sn,pn), M, "left" );
    cout << "Wiener predictor pn:   " << pn << endl;
    cout << "original\tpredicted\terror" << endl;
    Type realValue = sin( Type(M*TWOPI/10) );
    cout << setiosflags(ios::fixed) << setprecision(4)
         << realValue << "\t\t" << snPred[M-1] << "\t\t"
         << abs(realValue-snPred[M-1]) << endl << endl;

    return 0;
}
