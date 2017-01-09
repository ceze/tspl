/*****************************************************************************
 *                               fftmr_test.cpp
 *
 * Mixed Radix Algorithm FFT testing.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#include <iostream>
#include <iomanip>
#include <fftmr.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     LENGTH = 32;


int main()
{
    int     i, j, index, rows = LENGTH/4;

    Vector<Type> xn(LENGTH);
    Vector< complex<Type> > yn(LENGTH),
                              Xk(LENGTH);
    FFTMR<Type> Fourier;

    cout << "The original signal is: " << endl;
    for( i=0; i<rows; i++ )
    {
        cout << endl;
        for( j=0; j<3; j++ )
        {
            index = 3*i+j;
            xn[index] = i+j;
            cout << setiosflags(ios::fixed) << setprecision(6);
            cout << "\t" << xn[index];
        }
    }
    cout << endl << endl;

    Fourier.fft( xn, Xk );

    cout << "The Fourier transform of original signal is:" << endl;
    for( i=0; i<rows; i++ )
    {
        cout << endl;
        for( j=0; j<3; j++ )
        {
            index = 3*i+j;
            cout << setiosflags(ios::fixed) << setprecision(6);
            cout << "\t" << Xk[index];
        }
    }
    cout << endl << endl;

    Fourier.ifft( Xk, xn );
    cout << "The inverse Fourier transform is" << endl;
    for( i=0; i<rows; i++ )
    {
        cout << endl;
        for( j=0; j<3; j++ )
        {
            index = 3*i+j;
            cout << setiosflags(ios::fixed) << setprecision(6);
            cout << "\t" << xn[index];
        }
    }
    cout << endl << endl;

    cout << "The original signal is: " << endl;
    for( i=0; i<rows; i++ )
    {
        cout << endl;
        for( j=0; j<3; j++ )
        {
            index = 3*i+j;
            yn[index] = complex<double>(i,j);
            cout << setiosflags(ios::fixed) << setprecision(6);
            cout << "\t" << yn[index];
        }
    }
    cout << endl << endl;

    Fourier.fft( yn );
    cout << "The Fourier transform of original signal is:" << endl;
    for( i=0; i<rows; i++ )
    {
        cout << endl;
        for( j=0; j<3; j++ )
        {
            index = 3*i+j;
            cout << setiosflags(ios::fixed) << setprecision(6);
            cout << "\t" << yn[index];
        }
    }
    cout << endl << endl;

    Fourier.ifft( yn );
    cout << "The inverse Fourier transform is" << endl;
    for( i=0; i<rows; i++ )
    {
        cout << endl;
        for( j=0; j<3; j++ )
        {
            index = 3*i+j;
            cout << setiosflags(ios::fixed) << setprecision(6);
            cout << "\t" << yn[index];
        }
    }

    cout << endl << endl;
    return 0;
}
