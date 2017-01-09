/*****************************************************************************
 *                               lsfit_test.cpp
 *
 * Least square fitting testing.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <lsfitting.h>


using namespace std;
using namespace splab;


typedef double  Type;


int main()
{
    Type    A = 4.0,
            alpha = 1.0,
            beta = -2.0,
            gamma = -4.0,
            tmp = 0.0;

    int M = 100;
    Vector<Type> x = linspace( 0.01, 0.5*PI, M );
    Vector<Type> y(M);
    for( int i=0; i<M; ++i )
    {
        tmp = A * pow(x[i],alpha) * exp(beta*x[i]+gamma*x[i]*x[i]);
        y[i] = log(max(tmp,EPS));
    }

    Funcs<Type> phi;
    LSFitting<Type> lsf( x, y, phi );
    lsf.calcCoefs();
    Vector<Type> parms = lsf.getCoefs();
    parms(1) = exp(parms(1));

    cout << "The original parameters are:" << endl
         << A << endl << alpha << endl << beta << endl << gamma << endl << endl;
    cout << "The fitted parameters are:   "  << parms << endl;

    return 0;
}
