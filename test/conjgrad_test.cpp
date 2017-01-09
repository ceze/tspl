/*****************************************************************************
 *                               conjgrad_test.cpp
 *
 * Conjugate gradient optimal method testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <objfunc.h>
#include <conjgrad.h>


using namespace std;
using namespace splab;


typedef double  Type;


int main()
{
    Type a = 1.0,
         b = -1.0,
         c = -1.0;
    ObjFunc<Type> f( a, b, c );
    Vector<Type> x0(2);
    x0(1) = 0.5;
    x0(2) = 0.5;

    Type tolErr = 1.0e-6;
    ConjGrad< Type, ObjFunc<Type> > prp;
    prp.optimize( f, x0, x0.dim(), tolErr );
    if( prp.isSuccess() )
    {
        Vector<Type> xmin = prp.getOptValue();
        int N = prp.getItrNum();
        cout << "The iterative number is:   " << N << endl << endl;
        cout << "The number of function calculation is:   "
             << prp.getFuncNum() << endl << endl;
        cout << setiosflags(ios::fixed) << setprecision(4);
        cout << "The optimal value of x is:   " << xmin << endl;
        cout << "The minimum value of f(x) is:   " << f(xmin) << endl << endl;
        cout << "The gradient's norm at x is:   "
             << prp.getGradNorm()[N] << endl << endl;
    }
    else
        cout << "The optimal solution  cann't be found!" << endl;

    return 0;
}
