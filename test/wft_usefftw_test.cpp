/*****************************************************************************
 *                               wft_usefftw_test.cpp
 *
 * Windowed Fourier transform by using FFTW testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <wft_usefftw.h>
#include <vectormath.h>
#include <timing.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     Lg = 128;
const   int     Ls = 1000;
const   int     Fs = 1000;


int main()
{
	/******************************* [ signal ] ******************************/
	Type a = 0;
	Type b = Ls-1;
	Vector<Type> t = linspace(a,b,Ls) / Type(Fs);
	Vector<Type> s = sin( Type(400*PI) * pow(t,Type(2.0)) );

	/******************************** [ widow ] ******************************/
	a = 0;
	b = Type(Lg-1);
	Type u = (Lg-1)/Type(2);
	Type r = Lg/Type(8);
	t = linspace(a,b,Lg);
	Vector<Type> g = gauss(t,u,r);
	g = g/norm(g);

	/********************************* [ WFT ] *******************************/
	Type runtime = 0;
	Timing cnt;
	cout << "Taking windowed Fourier transform." << endl;
	cnt.start();
    Matrix< complex<Type> > coefs = wftFFTW( s, g );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	/******************************** [ IWFT ] *******************************/
	cout << "Taking inverse windowed Fourier transform." << endl;
	cnt.start();
	Vector<Type> x = iwftFFTW( coefs, g );
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	cout << "The relative error is : " << "norm(s-x) / norm(s) = "
		 << norm(s-x)/norm(s) << endl << endl;

	return 0;
}
