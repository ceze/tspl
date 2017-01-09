/*****************************************************************************
 *                               statistics_test.cpp
 *
 * Statistics routines testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <random.h>
#include <statistics.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 1000;


int main()
{
    Vector<Type> xn = randn( 37, Type(0.0), Type(1.0), N );

    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "The minimum maximum and median value of the sequence." << endl;
    cout << min(xn) << endl << max(xn) << endl << mid(xn) <<endl << endl;
    cout << "The mean, variance and standard variance of the sequence." << endl;
    cout << mean(xn) << endl << var(xn) << endl << stdVar(xn) << endl << endl;

    Vector<Type> yn = xn - mean(xn);
    cout << "The skew and kurtosis of the sequence." << endl;
    cout << skew(yn) << endl << kurt(yn) << endl << endl;

    cout << "The PDF of the sequence." << endl;
    cout << pdf(yn) << endl;

    return 0;
}
