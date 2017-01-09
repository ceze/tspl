/*****************************************************************************
 *                               window_test.cpp
 *
 * Windows function testing.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <window.h>


using namespace std;
using namespace splab;


typedef double  Type;


int main()
{
    int N = 5;
    Type A = 1.0;

    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "Rectangle window : " << rectangle(N,A) << endl;
    cout << "Bartlett window : " << bartlett(N,A) << endl;
    cout << "Hanning window : " << hanning(N,A) << endl;
    cout << "Hamming window : " << hamming(N,A) << endl;
    cout << "Blackman window : " << blackman(N,A) << endl;
    cout << "Kaiser window : " << kaiser(N,Type(8/PI),A) << endl;
    cout << "Gauss window : " << gauss(N,Type(2.5),A) << endl;

    cout << "Rectangle window : " << window("Rectangle",N,A) << endl;
    cout << "Bartlett window : " << window("Bartlett",N,A) << endl;
    cout << "Hanning window : " << window("Hanning",N,A) << endl;
    cout << "Hamming window : " << window("Hamming",N,A) << endl;
    cout << "Blackman window : " << window("Blackman",N,A) << endl;
    cout << "Kaiser window : " << window("Kaiser",N, Type(8/PI),A) << endl;
    cout << "Gauss window : " << window("Gauss",N,Type(2.5),A) << endl;

    return 0;
}
