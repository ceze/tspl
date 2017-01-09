/*****************************************************************************
 *                               time_test.cpp
 *
 * Timing function testing.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <cmath>
#include <timing.h>


using namespace std;
using namespace splab;


int main()
{
    double  var = 1;
    double  runTime = 0;
    Timing  time;

    time.start();
    for( int i=0; i<10000; ++i )
        for( int j=0; j<10000; ++j )
            var = sqrt(1.0*i*j);
    time.stop();
    runTime = time.read();

    cout << "The running time is : " << runTime << endl << endl;

    return 0;
}
