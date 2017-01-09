/*****************************************************************************
 *                               cwt_test.cpp
 *
 * Continuous wavelet transform testing.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <vectormath.h>
#include <statistics.h>
#include <timing.h>
#include <cwt.h>


using namespace std;
using namespace splab;


const   int     Ls = 1000;
const   double  fs = 1000.0;


int main()
{

	/******************************* [ signal ] ******************************/
	Vector<double> t = linspace( 0.0, (Ls-1)/fs, Ls );
	Vector<double> st = sin( 200*PI*pow(t,2.0) );
	st = st-mean(st);

	/******************************** [ CWT ] ********************************/
	Matrix< complex<double> > coefs;
	CWT<double> wavelet("morlet");
	wavelet.setScales( fs, fs/Ls, fs/2 );
	Timing cnt;
	double runtime = 0.0;
	cout << "Taking continuous wavelet transform(Morlet)." << endl;
	cnt.start();
	coefs = wavelet.cwtC(st);
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	/******************************** [ ICWT ] *******************************/
	cout << "Taking inverse continuous wavelet transform." << endl;
	cnt.start();
	Vector<double> xt = wavelet.icwtC(coefs);
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	cout << "The relative error is : " << endl;
	cout << "norm(st-xt) / norm(st) = " << norm(st-xt)/norm(st) << endl;
	cout << endl << endl;


    /******************************* [ signal ] ******************************/
    Vector<float> tf = linspace( float(0.0), (Ls-1)/float(fs), Ls );
	Vector<float> stf = sin( float(200*PI) * pow(tf,float(2.0) ) );
	stf = stf-mean(stf);

	/******************************** [ CWT ] ********************************/
	CWT<float> waveletf("mexiHat");
	waveletf.setScales( float(fs), float(fs/Ls), float(fs/2), float(0.25) );
	runtime = 0.0;
	cout << "Taking continuous wavelet transform(Mexican Hat)." << endl;
	cnt.start();
	Matrix<float> coefsf = waveletf.cwtR(stf);
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	/******************************** [ ICWT ] *******************************/
	cout << "Taking inverse continuous wavelet transform." << endl;
	cnt.start();
	Vector<float> xtf = waveletf.icwtR(coefsf);
	cnt.stop();
	runtime = cnt.read();
	cout << "The running time = " << runtime << " (ms)" << endl << endl;

	cout << "The relative error is : " << endl;
	cout << "norm(st-xt) / norm(st) = " << norm(stf-xtf)/norm(stf) << endl << endl;

	return 0;
}
