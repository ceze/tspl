/*****************************************************************************
 *                               utilities_test.cpp
 *
 * Utilities testing.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <string>
#include <utilities.h>


using namespace std;
using namespace splab;


const int N = 5;


int main()
{

	Vector<int> v1(N);
	for( int i=1; i<=v1.dim(); i++ )
		v1(i) = i;
	cout << "vector v1 : " << v1 << endl;
	Vector<int> v2(N);
	for( int i=1; i<=v2.dim(); ++i )
		v2(i) = i+N;
	cout << "vector v2 : " << v2 << endl;

	int N = 11;
	double a = 0;
	double b = 1.0;
	Vector<double> x = linspace( a, b, N );
	cout << N << " points linearly spaced from 0 to 1.0"
	     << x << endl;

	cout << "Flipping vector v1 from left to right : " << flip(v1) << endl;
	cout << "Shift vector v1 from left to right : " << shift(v1,2) << endl;
	cout << "Shift vector v1 from right to left : " << shift(v1,-2) << endl;
	cout << "Circle shift vector v1 from left to right : "
         << circshift(v1,2) << endl;
	cout << "Circle shift vector v1 from right to left : "
         << circshift(v1,-2) << endl;
	cout << "FFT shift of vector : " << fftshift(v1) << endl;

	cout << "Dyadic upsampling of vector v1 by zeros at the even position : "
		 << dyadUp( v1,0 ) << endl;
	cout << "Dyadic upsampling of vector v1 by zeros at the odd position : "
		 << dyadUp( v1,1 ) << endl;
	cout << "Dyadic downsampling of vector v1 by zeros at the even position : "
		 << dyadDown( v1,0 ) << endl;
	cout << "Dyadic downsampling of vector v1 by zeros at the odd position : "
		 << dyadDown( v1,1 ) << endl;

    Vector<float> sn(N);
    Vector< complex<float> > cn(N);
    for( int i=0; i<N; ++i )
    {
        sn[i] = float(sin(i*TWOPI/N));
        cn[i] = complex<float>( float(sin(i*TWOPI/N)), float(cos(i*TWOPI/N)) );
    }
    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "real signal sn : " << sn << endl;
    cout << "FFT interpolation of sn by factor fo 2 : "
		 << fftInterp( sn, 2 ) << endl;
    cout << "complex signal cn : " << cn << endl;
    cout << "FFT interpolation of sn by factor fo 2 : "
		 << fftInterp( cn, 2 ) << endl;
    cout << resetiosflags(ios::fixed);

	int n = 2;
	string dire = "left";
	string mode = "zpd";
	cout << "Extending vector v1 in left direction by zeros padding : "
		 << wextend( v1,n,dire,mode ) << endl;
	mode = "ppd";
	cout << "Extending vector v1 in left direction by periodic mode : "
		 << wextend( v1,n,dire,mode ) << endl;
	mode = "sym";
	cout << "Extending vector v1 in left direction by symmetric mode : "
		 << wextend( v1,n,dire,mode ) << endl;

	dire = "right";
	mode = "zpd";
	cout << "Extending vector v1 in right direction by zeros padding : "
		 << wextend( v1,n,dire,mode ) << endl;
	mode = "ppd";
	cout << "Extending vector v1 in right direction by periodic mode : "
		 << wextend( v1,n,dire,mode ) << endl;
	mode = "sym";
	cout << "Extending vector v1 in right direction by symmetric mode : "
		 << wextend( v1,n,dire,mode ) << endl;

	dire = "both";
	mode = "zpd";
	cout << "Extending vector v1 in both direction by zeros padding : "
		 << wextend( v1,n,dire,mode ) << endl;
	mode = "ppd";
	cout << "Extending vector v1 in both direction by periodic mode : "
		 << wextend( v1,n,dire,mode ) << endl;
	mode = "sym";
	cout << "Extending vector v1 in both direction by symmetric mode : "
		 << wextend( v1,n,dire,mode ) << endl;

	cout << "Keeping the center part of vector v1 : "
         << wkeep( v1,3,"center" ) << endl;
	cout << "Keeping the left part of vector v1 : "
         << wkeep( v1,3,"left" ) << endl;
	cout << "Keeping the right part of vector v1 : "
         << wkeep( v1,3,"right" ) << endl;
	cout << "Keeping the first(2) to first + L(3) elements of vector v1 : "
		 << wkeep( v1,3,2 ) << endl;

	cout << "The modulus of 2 divided by 5 is " << mod(2,5) << "." << endl;
	cout << "The modulus of -1 divided by 5 is " << mod(-1,5) << "." << endl;
	cout << endl;
	cout << "The nearest integer >= 10/2 is " << ceil(10,2) << "." << endl;
	cout << "The nearest integer >= 10/3 is " << ceil(10,3) << "." << endl;

	cout << endl;
	cout << "The numbers can be represented by the integer power of 2 "
		 << "from 0 to 1000 are : " << endl;

	return 0;
}
