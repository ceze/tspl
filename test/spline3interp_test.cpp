/*****************************************************************************
 *                               spline3interp_test.cpp
 *
 * Spline3 interpolation testing.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <spline3interp.h>


using namespace std;
using namespace splab;


typedef double   Type;


int main()
{
    // f(x) = 1 / (1+25*x^2)   -1 <= x <= 1
    Type xi[11] = { -1.0, -0.8, -0.6, -0.4, -0.2,
                    0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
    Type yi[11] = { 0.0385, 0.0588, 0.1, 0.2, 0.5,
                    1.0, 0.5, 0.2, 0.1, 0.0588, 0.0385 };
    Type Ml = 0.2105, Mr = 0.2105;
    Vector<Type> x( 11, xi ), y( 11, yi );

    Spline3Interp<Type> poly( x, y, Ml, Mr );
    poly.calcCoefs();

    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "Coefficients of cubic spline interpolated polynomial:   "
         << poly.getCoefs() << endl;

    cout << "The true and interpolated values:" << endl;
    cout << "0.0755" << "   " << poly.evaluate(-0.7) << endl
         << "0.3077" << "   " << poly.evaluate(0.3) << endl
         << "0.0471" << "   " << poly.evaluate(0.9) << endl << endl;

    return 0;
}
