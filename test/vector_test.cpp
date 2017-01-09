/*****************************************************************************
 *                               vector_test.cpp
 *
 * Vector class testing.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector.h>


using namespace std;
using namespace splab;


const int M = 3;


void display( const int *p, int length )
{
	for( int i=0; i<length; ++i )
		cout << p[i] << "\t" ;
	cout << endl;
}


int main()
{
	int k;
	int arrays[3] = {1,2,3};

	Vector<int> v1( M,arrays );
	k = 1;
	Vector<int> v2( M,k );

	Vector<int> v3( M );
	k = 0;
	v3 = k;

	Vector<int> v4( v1 );

	cout << "vector v1 : " << v1 << endl;
	cout << "vector v2 : " << v2 << endl;
	cout << "vector v3 : " << v3 << endl;
	cout << "vector v4 : " << v4 << endl;

	display( v4, M );
	cout << endl;

	v4.resize( 5 );
	Vector<int>::iterator itr = v4.begin();
	while( itr != v4.end() )
        *itr++ = 1;
	cout << "new vector v4 : " << v4 << endl;
	v4 = v1;
	cout << "new vector v4 : " << v4 << endl;

	k = 2;
	v3 = k+v1;
	cout << "v3 = k + v1 : " << v3 << endl;
	v3 += k;
	cout << "v3 += k : " << v3 << endl;

    v3 = v1-k;
	cout << "v3 = v1 - k : " << v3 << endl;
	v3 = k-v1;
	cout << "v3 = k - v1 : " << v3 << endl;
	v3 -= k;
	cout << "v3 -= k : " << v3 << endl;

	v3 = k*v1;
	cout << "v3 = k * v1 : " << v3 << endl;
	v3 *= k;
	cout << "v3 *= k : " << v3 << endl;

	v3 = v1/k;
	cout << "v3 = v1 / k : " << v3 << endl;
	v3 = k/v1;
	cout << "v3 = k / v1 : " << v3 << endl;
	v3 /= k;
	cout << "v3 /= k : " << v3 << endl;

    v3 = v1+v2;
	cout << "v3 = v1 + v2 : " << v3 << endl;
	v3 += v1;
	cout << "v3 += v1 : " << v3 << endl;

	v3 = v1-v2;
	cout << "v3 = v1 - v2 : " << v3 << endl;
	v3 -= v1;
	cout << "v3 -= v1 : " << v3 << endl;

	v3 = v1*v2;
	cout << "v3 = v1 * v2 : " << v3 << endl;
	v3 *= v1;
	cout << "v3 *= v1 : " << v3 << endl;

	v3 = v1/v2;
	cout << "v3 = v1 / v2 : " << v3 << endl;
	v3 /= v1;
	cout << "v3 /= v1 : " << v3 << endl;

    cout << "minimum element of v1 :   " << min(v1) << endl << endl;
    cout << "maximum element of v1 :   " << max(v1) << endl << endl;
    cout << "L2 norm of v3 :   " << norm( v3 ) << endl << endl;
	cout << "inner product of v1 and v2 :   "
         << dotProd( v1, v2 ) << endl << endl;

	complex<double> z = -1.0;
	Vector< complex<double> > v( M );
	v[0] = polar( 1.0,PI/4 );
	v[1] = polar( 1.0,PI );
	v[2] = complex<double>( 1.0,-1.0 );
	Vector< complex<double> > u = v*z;

    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "convert from real to complex vector: "
         << complexVector(v3) << endl;
    cout << "convert from real to complex vector: "
         << complexVector(v3,-v1) << endl;
	cout << "complex vector v : " << v << endl;
	cout << "complex vector u = -v : " << u << endl;
	cout << "norm of coplex vector v : " << norm(v) << endl << endl;
	cout << "dot product of complex vector v and 1+u: "
         << dotProd(v,u-z) << endl << endl;

    int N = 5;
    Vector<double> x = linspace( 0.0, TWOPI, N );
    Vector< complex<float> > cv(N);
    for( int i=0; i<N; ++i )
        cv[i] = complex<float>( float(sin(x[i])), float(cos(x[i])) );
    cout << "Complex vector vc : " << cv << endl;
    cout << "Absolute of vc : " << abs(cv) << endl;
    cout << "Angle of vc : " << arg(cv) << endl;
    cout << "Real part of vc : " << real(cv) << endl;
    cout << "Imaginary part of vc : " << imag(cv) << endl;

    cout << resetiosflags(ios::fixed);
	Vector< Vector<double> > v2d1( M );
	for( int i=0; i<M; ++i )
	{
		v2d1[i].resize( M );
		for( int j=0; j<M; ++j )
			v2d1[i][j] = double( i+j );
	}

	cout << "two dimension vector v2d1 : " << endl;
	Vector< Vector<double> >::const_iterator itrD2 = v2d1.begin();
	int rowNum = 0;
	while( itrD2 != v2d1.end() )
	    cout << "the " << rowNum++ << "th row : " << *itrD2++ << endl;

	Vector< Vector<double> > v2d2 = v2d1+v2d1;
	cout << "two dimension vector v2d2 = v2d1 + v2d1 : " << v2d2;


	return 0;
}
