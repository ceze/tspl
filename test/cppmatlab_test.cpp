/*****************************************************************************
 *                             CppMatlab_test.cpp
 *
 * C++ and Matlab mixed programming testing.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <cstring>
#include <vectormath.h>
#include <matrixmath.h>
#include <wft.h>
#include "engine.h"


using namespace std;
using namespace splab;


typedef double  Type;
const   int     Lg = 128;
const   int     Ls = 1000;
const   Type    Fs = 1000;


int main()
{
    /******************************* [ signal ] ******************************/
	Vector<Type> t = linspace( Type(0), Type(Ls-1), Ls ) / Type(Fs);
	Vector<Type> s = sin( Type(400*PI) * pow(t,Type(2.0)) );

	/******************************** [ widow ] ******************************/
	t = linspace(Type(0),Type(Lg-1),Lg);
	Vector<Type> g = gauss( t, (Lg-1)/Type(2), Lg/Type(8) );

	/********************************* [ WFT ] *******************************/
	cout << "Taking windowed Fourier transform." << endl;
    Matrix< complex<Type> > coefs = wft( s, g );

	/******************************** [ IWFT ] *******************************/
	cout << "Taking inverse windowed Fourier transform." << endl;
	Vector<Type> x = iwft( coefs, g );

	cout << "The relative error is : " << "norm(s-x) / norm(s) = "
		 << norm(s-x)/norm(s) << endl << endl;

    /******************************** [ PLOT ] *******************************/
	Engine *ep  = engOpen( NULL );
	if( !ep )
	{
		cerr << "Cannot open Matlab Engine!" << endl;
		exit(1);
	}

    Matrix<Type> C = trT( abs(coefs) );
    int M = C.rows(), N = C.cols();

    // define mxArray as 1-by-1 Real Scalar
	mxArray *mFs = mxCreateDoubleMatrix( 1, 1, mxREAL );

    // define mxArray as N-by-1 Real Vector
	mxArray *ms = mxCreateDoubleMatrix( Ls, 1, mxREAL );
    mxArray *mx = mxCreateDoubleMatrix( Ls, 1, mxREAL );

	// define mxArray as N-by-M Real Matrix, BECAUSE matalb is ROW MAJOR
	// and C/C++ is COLUMN MAJOR, the row of SP++ matrix is copied as the
	// column of Matlab matrix
	mxArray *mC = mxCreateDoubleMatrix( N, M, mxREAL );

    // array copy from Scalar to mxArray.
	memcpy( mxGetPr(mFs), &Fs, sizeof(Type) );

    // array copy from Vectors to mxArray.
	memcpy( mxGetPr(ms), s, Ls*sizeof(Type) );
	memcpy( mxGetPr(mx), x, Ls*sizeof(Type) );

	// array copy from Matrix to mxArray.
	memcpy( mxGetPr(mC), C, C.size()*sizeof(Type) );

    // send command to Matlab engine
    engPutVariable( ep, "fs", mFs );
	engPutVariable( ep, "s", ms );
	engPutVariable( ep, "x", mx );
	engPutVariable( ep, "C", mC );

    // Matlab commands
	const char *mCmd =  " \
        figure('name','C++ and Matlab Mixed Programming Testing'); \
        hFN = floor(size(C,1)/2);   tN = size(C,2); \
        subplot(2,2,1); plot((0:tN-1), s); \
        axis([0,tN,min(s),max(s)]); \
        xlabel('Time (ms)', 'fontsize',12); ylabel('Amplitude', 'fontsize',12); \
        title('(a)'); \
        subplot(2,2,2); pcolor((0:tN-1),(0:hFN)'/hFN, C(1:hFN+1,:)); \
        shading interp; \
        axis([0,tN, 0,1]); \
        yt = 0 : 0.2 : 1; set(gca, 'YTick',yt); set(gca, 'YTickLabel',fs*yt/2); \
        xlabel('Time (ms)','fontsize',12); ylabel('Frequency (Hz)','fontsize',12); \
        title('(b)'); \
        subplot(2,2,3); plot((0:tN-1),x); \
        axis([0,tN,min(x),max(x)]); \
        xlabel('Time (ms)','fontsize',12); \
        ylabel('Amplitude', 'fontsize',12); \
        title('(c)'); \
        subplot(2,2,4); e=s-x; plot((0:tN-1),e); \
        axis([0,tN,min(e),max(e)]); \
        xlabel('Time (ms)','fontsize',12); \
        ylabel('Amplitude', 'fontsize',12); \
        title('(d)'); \
        ";

    // send command to Matlab engine
    engEvalString( ep, mCmd );

	// delete mxArray
	mxDestroyArray( mFs );
	mxDestroyArray( ms );
	mxDestroyArray( mx );
	mxDestroyArray( mC );
    system( "pause" );
	engClose(ep);

	return 0;
}
