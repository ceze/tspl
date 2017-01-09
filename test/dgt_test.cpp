/*****************************************************************************
 *                               dgt_test.cpp
 *
 * Discrete Gabor transform testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <vectormath.h>
#include <timing.h>
#include <dgt.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     Fs = 1000;
const   int     Ls = 10000;
const   int     Lg = 80;
const   int     N  = 40;
const   int     dM = 10;      // over sampling ratio is N/dM


int main()
{

	/******************************* [ signal ] ******************************/
	Type a = 0;
	Type b = Ls-1;
	Vector<Type> t = linspace( a, b, Ls ) / Type(Fs);
	Vector<Type> st = cos( Type(400*PI) * pow(t,Type(2.0)) );

	/******************************* [ widow ] ******************************/
	a = 0;
	b = Type( Lg-1 );
	Type r = sqrt( dM*N / Type(2*PI) );
	Type u = (Lg-1) / Type(2.0);
	t = linspace( a, b, Lg );
	Vector<Type> h = gauss( t, u, r );
	h = h / norm(h);

	/***************************** [ daul function ] *************************/
	Type runtime = 0.0;
	Timing cnt;
	cout << "Compute daul function." << endl;
	cnt.start();
	Vector<Type> g = daul( h, N, dM );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	/******************************** [ DGT ] ********************************/
	cout << "Taking discrete Gabor transform." << endl;
	cnt.start();
	Matrix< complex<Type> > C = dgt( st, g, N, dM, "sym" );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	/******************************** [ IDGT ] *******************************/
	cout << "Taking inverse discrete Gabor transform." << endl;
	cnt.start();
	Vector<Type> xt = idgt( C, h, N, dM );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	cout <<"The relative error : norm(s-x) / norm(s) = "
	     << norm(st-xt)/norm(st) << endl << endl;

	return 0;
}
