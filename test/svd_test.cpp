/*****************************************************************************
 *                               svd_test.cpp
 *
 * SVD class testing.
 *
 * Zhang Ming, 2010-01 (revised 2010-08), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <svd.h>


using namespace std;
using namespace splab;


typedef double  Type;


int main()
{
//    Matrix<Type> A(4,4);
//	A(1,1) = 16;    A(1,2) = 2;     A(1,3) = 3;     A(1,4) = 13;
//	A(2,1) = 5;     A(2,2) = 11;    A(2,3) = 10;    A(2,4) = 8;
//	A(3,1) = 9;     A(3,2) = 7;     A(3,3) = 6;     A(3,4) = 12;
//	A(4,1) = 4;     A(4,2) = 14;    A(4,3) = 15;    A(4,4) = 1;

//    Matrix<Type> A(4,2);
//	A(1,1) = 1;     A(1,2) = 2;
//	A(2,1) = 3;     A(2,2) = 4;
//	A(3,1) = 5;     A(3,2) = 6;
//	A(4,1) = 7;     A(4,2) = 8;

    Matrix<Type> A(2,4);
	A(1,1) = 1;     A(1,2) = 3;     A(1,3) = 5;     A(1,4) = 7;
	A(2,1) = 2;     A(2,2) = 4;     A(2,3) = 6;     A(2,4) = 8;

	SVD<Type> svd;
	svd.dec(A);

	Matrix<Type> U = svd.getU();
	Matrix<Type> V = svd.getV();
	Matrix<Type> S = svd.getSM();

    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "Matrix--A: " << A << endl;
	cout << "Matrix--U: " << U << endl;
	cout << "Vector--S: " << S << endl;
	cout << "Matrix--V: " << V << endl;
    cout << "Matrix--A - U * S * V^T:  "
         << A- U*S*trT(V) << endl;
//         << A- U*multTr(S,V) << endl;

	cout << "The rank of A : " << svd.rank() << endl << endl;
	cout << "The condition number of A : " << svd.cond() << endl << endl;

	return 0;
}
