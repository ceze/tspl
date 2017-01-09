/*****************************************************************************
 *                               nle_test.cpp
 *
 * Rooting of nonlinear equations testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <iomanip>
#include <nleroots.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 2;


int main()
{
    Vector<Type> X0(N);
    cout << setiosflags(ios::fixed) << setprecision(8);

    NLEqus<Type> G;
    X0(1) = 0;  X0(2) = 1;
    cout << "Seidel iteration method :   " << seidel( G, X0 ) << endl;

    NLFuncs<Type> F;
	X0(1) = 2;  X0(2) = 0;
	cout << "Newton iteration method :   " << newton( F, X0 ) << endl;

	return 0;
}
