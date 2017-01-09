/*****************************************************************************
 *                             vectormath_test.cpp
 *
 * Math functions of vector testing.
 *
 * Zhang Ming, 2010-08, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <vectormath.h>


using namespace std;
using namespace splab;


int main()
{

	int     N = 5;
	double  a = 0,
            b = 2*PI;
	Vector<double> x = linspace( a, b, N );

	cout << setiosflags(ios::fixed) << setprecision(4);
	cout << x << endl;
	cout << "sin of x : " << sin(x) << endl;
	cout << "cos of x : "<< cos(x) << endl;
	cout << "tan of x : " << tan(x) << endl;
	cout << "asin of x : "<< asin(x) << endl;
	cout << "acos of x : " << acos(x) << endl;
	cout << "atan of x : " << atan(x) << endl;

	cout << "exp of x : "<< exp(x) << endl;
	cout << "log of x : " << log(x) << endl;
	cout << "log10 of x : " << log10(x) << endl;

    a = 2.0;
	cout << "sqrt of x : "<< sqrt(x) << endl;
	cout << "pow of x : " << pow(x,x) << endl;
	cout << "pow of x : " << pow(x,a) << endl;
	cout << "pow of x : " << pow(a,x) << endl;

	cout << "The standard normal distribution : " << gauss( x, a, b ) << endl;

	return 0;
}
