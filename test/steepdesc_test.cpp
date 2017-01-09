/*****************************************************************************
 *                               steepdesc_test.cpp
 *
 * Steepest descent method testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <objfunc.h>
#include <steepdesc.h>


using namespace std;
using namespace splab;


typedef float  Type;


int main()
{
    Type a = 1.0,
         b = -1.0,
         c = -1.0;
    ObjFunc<Type> f( a, b, c );
    Vector<Type> x0(2);
    x0(1) = Type(0.0);
    x0(2) = Type(0.0);

    Type tolErr = 1.0e-3;
    SteepDesc< Type, ObjFunc<Type> > steep;
    steep.optimize( f, x0, tolErr );
    if( steep.isSuccess() )
    {
        Vector<Type> xmin = steep.getOptValue();
        int N = steep.getItrNum();
        cout << "The iterative number is:   " << N << endl << endl;
        cout << "The number of function calculation is:   "
             << steep.getFuncNum() << endl << endl;
        cout << setiosflags(ios::fixed) << setprecision(4);
        cout << "The optimal value of x is:   " << xmin << endl;
        cout << "The minimum value of f(x) is:   " << f(xmin) << endl << endl;
        cout << "The gradient's norm at x is:   "
             << steep.getGradNorm()[N] << endl << endl;
    }
    else
        cout << "The optimal solution  can't be found!" << endl;

    return 0;
}
