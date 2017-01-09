/*****************************************************************************
 *                               nle_test.cpp
 *
 * Rooting of nonlinear equation testing.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <iomanip>
#include <nleroot.h>


using namespace std;
using namespace splab;


typedef double  Type;


int main()
{
    Type x01, x02, x03;
    NLFunc<Type> f( 1.0, -3.0, 1.0 );

    cout << setiosflags(ios::fixed) << setprecision(8);
    x01 = bisection( f, 0.0, 2.0, 1.0e-6 );
	cout << "Bisection method :   " << x01 << endl <<endl;

	x02 = newton( f, 0.0, 1.0e-6, 100 );
	cout << "Newton method :      " << x02 << endl <<endl;

	x03 = secant( f, 1.0, 3.0, 1.0e-6, 100 );
	cout << "Secant method :      " << x03 << endl <<endl;

	return 0;
}
