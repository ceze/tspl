/*****************************************************************************
 *                            classicalpse_test.cpp
 *
 * Classical power spectrum estimation testing.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <cstring>
#include <vectormath.h>
#include <classicalpse.h>
#include "engine.h"


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 100;
const   int     K = 25;
const   int     M = 50;
const   int     L = 200;


int main()
{
    /******************************* [ signal ] ******************************/
    int mfn = L/2+1;
	Type amp1 = Type(1.0),
         amp2 = Type(0.5);
    Type f1 = Type(0.3),
         f2 = Type(0.4);
    Vector<Type> tn = linspace(Type(0), Type(N-1), N );
    Vector<Type> sn = amp1*sin(TWOPI*f1*tn) + amp2*sin(TWOPI*f2*tn);

	/******************************** [ widow ] ******************************/
	Vector<Type> wn = window( "Hamming", M, Type(1.0) );

	/********************************* [ PSD ] *******************************/
//	Vector<Type> Ps = correlogramPSE( sn, L );
//	Vector<Type> Ps = periodogramPSE( sn, wn, L );
//    Vector<Type> Ps = bartlettPSE( sn, M, L );
	Vector<Type> Ps = welchPSE( sn, wn, K, L );
//	Vector<Type> Ps = btPSE( sn, wn, L );

    /******************************** [ PLOT ] *******************************/
	Engine *ep  = engOpen( NULL );
	if( !ep )
	{
		cerr << "Cannot open Matlab Engine!" << endl;
		exit(1);
	}

	mxArray *msn = mxCreateDoubleMatrix( N, 1, mxREAL );
    mxArray *mPs = mxCreateDoubleMatrix( mfn, 1, mxREAL );
	memcpy( mxGetPr(msn), sn, N*sizeof(Type) );
	memcpy( mxGetPr(mPs), Ps, mfn*sizeof(Type) );
	engPutVariable( ep, "sn", msn );
	engPutVariable( ep, "Ps", mPs );

	const char *mCmd =  " figure('name','Welch Method of Spectrum Estimation'); \
        N = length(sn); mfn = length(Ps); \
        subplot(2,1,1); \
            plot((0:N-1), sn); \
            axis([0,N,min(sn),max(sn)]); \
            title('(a)   Signal', 'FontSize',12); \
            xlabel('Samples', 'FontSize',12); \
            ylabel('Amplitude', 'FontSize',12); \
        subplot(2,1,2); \
            h = stem((0:mfn-1)/(mfn-1)/2, Ps); \
            axis([0,0.5,min(Ps),max(Ps)]); \
            set(h,'MarkerFaceColor','blue'); \
            set(gca, 'XTick', 0:0.1:0.5); \
            grid on; \
            title('(b)   Spectrum', 'FontSize',12); \
            xlabel('Normalized Frequency ( f / fs )', 'FontSize',12); \
            ylabel('Amplitude', 'FontSize',12); ";
    engEvalString( ep, mCmd );

	mxDestroyArray( msn );
	mxDestroyArray( mPs );
    system( "pause" );
	engClose(ep);

	return 0;
}
