/*****************************************************************************
 *                            pseudoinverse_test.cpp
 *
 * Matrix pseudoinverse testing.
 *
 * Zhang Ming, 2010-08 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <pseudoinverse.h>


using namespace std;
using namespace splab;


typedef double  Type;


int main()
{
	Matrix<Type> A(8,6), invA;
    A[0][0]=64; A[0][1]=2;  A[0][2]=3;  A[0][3]=61; A[0][4]=60; A[0][5]=6;
    A[1][0]=9;  A[1][1]=55; A[1][2]=54; A[1][3]=12; A[1][4]=13; A[1][5]=51;
    A[2][0]=17; A[2][1]=47; A[2][2]=46; A[2][3]=20; A[2][4]=21; A[2][5]=43;
    A[3][0]=40; A[3][1]=26; A[3][2]=27; A[3][3]=37; A[3][4]=36; A[3][5]=30;
    A[4][0]=32; A[4][1]=34; A[4][2]=35; A[4][3]=29; A[4][4]=28; A[4][5]=38;
    A[5][0]=41; A[5][1]=23; A[5][2]=22; A[5][3]=44; A[5][4]=45; A[5][5]=19;
    A[6][0]=49; A[6][1]=15; A[6][2]=14; A[6][3]=52; A[6][4]=53; A[6][5]=11;
    A[7][0]=8;  A[7][1]=58; A[7][2]=59; A[7][3]=5;  A[7][4]=4;  A[7][5]=62;

    invA = pinv(A);
	cout << setiosflags(ios::fixed) << setprecision(4);
	cout << "The original matrix A is : " << A << endl;
    cout << "The pseudoinverse matrix of A is : " << invA << endl;
    cout << "The multiplication of A and its inverse is : " << A*invA << endl;

    Matrix<complex<Type> > cA = complexMatrix(A,-A);
    Matrix<complex<Type> > cPIA = pinv(cA);
    cout << "The original complex matrix A is : "
         << setprecision(0) << cA << endl;
    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "The real part of the pseudoinverse matrix of A is: "
         << real(cPIA) << endl;
    cout << "The real imaginary of the pseudoinverse matrix of A is: "
         << imag(cPIA) << endl;

	return 0;
}
