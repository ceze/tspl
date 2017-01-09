/*****************************************************************************
 *                               wvd_test.cpp
 *
 * Wigner-Ville distribution testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <timing.h>
#include <wvd.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     Ls = 100;


int main()
{
	/******************************* [ signal ] ******************************/
	Vector<Type> t = linspace( Type(-0.5), Type(0.5-1.0/Ls), Ls );
	Vector< complex<Type> > sn(Ls);
	for( int i=0; i<Ls; ++i )
	{
	    Type tip2 = t[i]*t[i],
             freq = Type(TWOPI*25);
        sn[i] = exp(-32*tip2) * polar(Type(1),freq*t[i]) * polar(Type(1),freq*tip2);
	}
	cout << setiosflags(ios::fixed) << setprecision(4);

	/********************************* [ WVD ] *******************************/
	Type runtime = 0;
	Timing cnt;
	cout << "Computing Wigner-Wille distribution." << endl << endl;
	cnt.start();
    Matrix<Type> coefs = wvd( sn );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

    Vector<Type> timeMarg = mean(coefs),
                 freqMarg = mean(trT(coefs));
	cout << "The time marginal condition is: " << timeMarg << endl;
	cout << "The frequency marginal condition is: " << freqMarg << endl;

	return 0;
}
