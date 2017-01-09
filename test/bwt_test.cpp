/*****************************************************************************
 *                               bwt_test.cpp
 *
 * Dyadic wavelet transform testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <timing.h>
#include <bwt.h>


using namespace std;
using namespace splab;


const int Ls = 100;


int main()
{

	/************************** [ signal ] *************************/
	Vector<double> s(Ls);
	for(int i=0; i<Ls; i++)
	{
		if(i<Ls/4)
			s[i] = 0.0;
		else if(i<2*Ls/4)
			s[i] = 1.0;
		else if(i<3*Ls/4)
			s[i] = 3.0;
		else
			s[i] = 0.0;
	}

	/*************************** [ BWT ] ***************************/
	int level = 3;
	Timing cnt;
	double runtime = 0.0;
	cout << "Taking dyadic wavelet transform." << endl;
	cnt.start();
	Vector< Vector<double> > coefs = bwt( s, level );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (s)" << endl << endl;

	/*************************** [ IBWT ] **************************/
	cout << "Taking inverse dyadic wavelet transform." << endl;
	cnt.start();
	Vector<double> x = ibwt( coefs, level );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (s)" << endl << endl;

	cout << "The relative error is : norm(s-x) / norm(s) = "
         << norm(s-x)/norm(s) << endl << endl;

	return 0;
}
