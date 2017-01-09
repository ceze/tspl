/*****************************************************************************
 *                            parametricpse_test.cpp
 *
 * parametric power spectrum estimation testing.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <cstring>
#include <random.h>
#include <vectormath.h>
#include <parametricpse.h>
#include "engine.h"


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 50;
const   int     L = 200;
const   int     yuleOrder = 4;
const   int     burgOrder = 4;
const   int     lplsOrder = 4;


int main()
{
    /******************************* [ signal ] ******************************/
    cout << setiosflags(ios::fixed) << setprecision(4);
    int mfn = L/2+1;
    Type amp1 = Type(1.0),
         amp2 = Type(1.0);
    Type f1 = Type(0.2),
         f2 = Type(0.4);
    Type sigma2, SNR;

    Vector<Type> tn = linspace(Type(0), Type(N-1), N );
    Vector<Type> sn = amp1*sin(TWOPI*f1*tn) + amp2*sin(TWOPI*f2*tn);
    Vector<Type> wn = randn( 37, Type(0.0), Type(0.1), N );
    Vector<Type> xn = sn + wn;
    SNR = 20*log10(norm(sn)/norm(wn));
    cout << "The SNR = " << SNR << endl << endl;

    /********************************* [ PSD ] *******************************/
//    Vector<Type> ak = yulewalkerPSE( xn, yuleOrder, sigma2 );
//    cout << "The estimated AR coefficients by Youle-Walker method are: "
//         << ak << endl;

//    Vector<Type> ak = burgPSE( xn, burgOrder, sigma2 );
//    cout << "The estimated AR coefficients by Burg method are: "
//         << ak << endl;

    Vector<Type> ak = fblplsPSE( xn, lplsOrder, sigma2 );
    cout << "The estimated AR coefficients by Youle-Walker method are: "
         << ak << endl;

    cout << "The estimated variance is:   " << sigma2 << endl << endl;

    Vector<Type> bk(1,Type(1.0));
    Vector<Type> Px = armaPSD( ak, bk, sigma2, L );

    /******************************** [ PLOT ] *******************************/
	Engine *ep  = engOpen( NULL );
	if( !ep )
	{
		cerr << "Cannot open Matlab Engine!" << endl;
		exit(1);
	}

	mxArray *mxn = mxCreateDoubleMatrix( N, 1, mxREAL );
    mxArray *mPx = mxCreateDoubleMatrix( mfn, 1, mxREAL );
	memcpy( mxGetPr(mxn), xn, N*sizeof(Type) );
	memcpy( mxGetPr(mPx), Px, mfn*sizeof(Type) );
	engPutVariable( ep, "xn", mxn );
	engPutVariable( ep, "Px", mPx );

	const char *mCmd =  " figure('name','FBLPLS Method of Spectrum Estimation'); \
        N = length(xn); mfn = length(Px); \
        subplot(2,1,1); \
            plot((0:N-1), xn); \
            axis([0,N,min(xn),max(xn)]); \
            title('(a)   Signal', 'FontSize',12); \
            xlabel('Samples', 'FontSize',12); \
            ylabel('Amplitude', 'FontSize',12); \
        subplot(2,1,2); \
            h = stem((0:mfn-1)/(mfn-1)/2, Px); \
            axis([0,0.5,min(Px),max(Px)]); \
            set(h,'MarkerFaceColor','blue'); \
            set(gca, 'XTick', 0:0.1:0.5); \
            grid on; \
            title('(b)   Spectrum', 'FontSize',12); \
            xlabel('Normalized Frequency ( f / fs )', 'FontSize',12); \
            ylabel('Amplitude', 'FontSize',12); ";
    engEvalString( ep, mCmd );

	mxDestroyArray( mxn );
	mxDestroyArray( mPx );
    system( "pause" );
	engClose(ep);

	return 0;
}
