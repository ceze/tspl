/*****************************************************************************
 *                               bfgs_test.cpp
 *
 * BFGS method testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <objfunc.h>
#include <bfgs.h>


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
    x0(1) = Type(0.0);
    x0(2) = Type(0.0);

    BFGS< Type, ObjFunc<Type> > bfgs;
    bfgs.optimize( f, x0 );
    if( bfgs.isSuccess() )
    {
        Vector<Type> xmin = bfgs.getOptValue();
        int N = bfgs.getItrNum();
        cout << "The iterative number is:   " << N << endl << endl;
        cout << "The number of function calculation is:   "
             << bfgs.getFuncNum() << endl << endl;
        cout << setiosflags(ios::fixed) << setprecision(4);
        cout << "The optimal value of x is:   " << xmin << endl;
        cout << "The minimum value of f(x) is:   " << f(xmin) << endl << endl;
        cout << "The gradient's norm at x is:   "
             << bfgs.getGradNorm()[N] << endl << endl;
    }
    else
        cout << "The optimal solution  cann't be found!" << endl;

    return 0;
}
