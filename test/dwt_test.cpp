/*****************************************************************************
 *                               dwt_test.cpp
 *
 * Discrete wavelet transform testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <dwt.h>
#include <timing.h>


using namespace std;
using namespace splab;


typedef float   Type;
const   int     Ls = 1000;


int main()
{

	/******************************* [ signal ] ******************************/
	Vector<Type> s(Ls);
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

	/******************************** [ DWT ] ********************************/
	int level = 3;
	DWT<Type> discreteWT("db4");
	Timing cnt;
	double runtime = 0.0;
	cout << "Taking discrete wavelet transform." << endl;
	cnt.start();
	Vector<Type> coefs = discreteWT.dwt( s, level );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (s)" << endl << endl;

	/******************************** [ IDWT ] *******************************/
	level = 0;
    cout << "Taking inverse discrete wavelet transform." << endl;
	cnt.start();
	Vector<Type> x = discreteWT.idwt( coefs, level );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (s)" << endl << endl;

	cout << "The relative error is : norm(s-x) / norm(s) = "
	     << norm(s-x)/norm(s) << endl << endl;

	return 0;
}
