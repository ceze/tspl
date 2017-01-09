/*****************************************************************************
 *                            eigenanalysispse_test.cpp
 *
 * Eigenanalysis spectrum estimation testing.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <cstring>
#include <random.h>
#include <vectormath.h>
#include <eigenanalysispse.h>
#include "engine.h"


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 200;
const   int     M = 20;
const   int     L = 200;


int main()
{
    /******************************* [ signal ] ******************************/
    cout << setiosflags(ios::fixed) << setprecision(4);
    int mfn = L/2+1;
    Type amp1 = Type(1.0),
         amp2 = Type(1.0);
    Type f1 = Type(0.2),
         f2 = Type(0.25);
    Type SNR;

    Vector<Type> tn = linspace(Type(0), Type(N-1), N );
    Vector<Type> sn = amp1*sin(Type(TWOPI)*f1*tn) + amp2*sin(Type(TWOPI)*f2*tn);
    Vector<Type> wn = randn( 37, Type(0.0), Type(0.2), N );
    Vector<Type> xn = sn + wn;
    SNR = 20*log10(norm(sn)/norm(wn));
    cout << "The SNR = " << SNR << endl << endl;

    /********************************* [ PSD ] *******************************/
    int p = orderEst( xn, M );
    p = p/2*2;
    cout << "The model order that minimizes the MDL criterion is:   p = "
         << p << endl << endl;

//    Vector<Type> Px = caponPSE( xn, M, L );

//    Vector<Type> Px = pisarenkoPSE( xn, M, p, L );

//    Vector<Type> Px = musicPSE( xn, M, p, L );

    Vector<Type> fk = espritPSE( xn, M, p );
    Vector<Type> Px(mfn);
    for( int k=0; k<p; ++k )
    {
        int index = int( L*abs(fk[k]) + 0.5 );
        if( index != 0 && index != L/2 )
        Px[index] = Type(1.0);
    }

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
            set(gca, 'XTick', 0:0.05:0.5); \
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
