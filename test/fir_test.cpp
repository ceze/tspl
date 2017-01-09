/*****************************************************************************
 *                               fir_test.cpp
 *
 * FIR class testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <fir.h>


using namespace std;
using namespace splab;


int main()
{
    string  wType = "Hamming";

//    string  fType = "lowpass";
//    double  fs = 1000,
//            fpass = 200,
//            apass = -3,
//            fstop = 300,
//            astop = -20;
//    FIR fir( fType, wType );
//    fir.setParams( fs, fpass, apass, fstop, astop );

//    string  fType = "highpass";
//    double  fs = 1000,
//            fstop = 200,
//            astop = -20,
//            fpass = 300,
//            apass = -3;
//    FIR fir( fType, wType );
//    fir.setParams( fs, fstop, astop, fpass, apass );

//    string  fType = "bandpass";
//    double  fs = 1000,
//            fstop1 = 100,
//            astop1 = -20,
//            fpass1 = 200,
//            fpass2 = 300,
//            apass1 = -3,
//            fstop2 = 400,
//            astop2 = -20;
//    FIR fir( fType, wType );
//    fir.setParams( fs, fstop1, astop1, fpass1, fpass2, apass1, fstop2, astop2 );

    string  fType = "bandstop";
    double  fs = 1000,
            fpass1 = 100,
            apass1 = -3,
            fstop1 = 200,
            fstop2 = 300,
            astop1 = -20,
            fpass2 = 400,
            apass2 = -3;
    FIR fir( fType, wType );
    fir.setParams( fs, fpass1, apass1, fstop1, fstop2, astop1, fpass2, apass2 );

    fir.design();
    fir.dispInfo();

    cout << endl;
    return 0;
}
