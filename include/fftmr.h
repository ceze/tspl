/*
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 or any later version.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details. A copy of the GNU General Public License is available at:
 * http://www.fsf.org/licensing/licenses
 */


/*****************************************************************************
 *                                fftmr.h
 *
 * Fast Fourier Transform with Mixed Radix Algorithm
 *
 * This class is designed for calculating discrete Fourier transform and
 * inverse discrete Fourier transform of 1D signals by using Mixed-Radix
 * Algorithms. The length of signals should equal to powers of 2.
 *
 * The algorithm is modified from "kube-gustavson-fft.c". Which is a pretty
 * fast FFT algorithm. Not super-optimized lightning-fast, but very good
 * considering its age and relative simplicity.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef FFTMR_H
#define FFTMR_H


#include <vector.h>


namespace splab
{

    /**
     * complex node converted from "C" version for this algorithm
     */
    template<typename Type>
    struct Complex
    {
        Type re;
        Type im;
    };


    /**
     * two routines frequently used in FFT algorithm
     */
    bool    isPower2( int );
    int     fastLog2( int );


    /**
     * Radix based FFT class
     */
	template<typename Type>
	class FFTMR
	{

	public:

		FFTMR();
		~FFTMR();

        void fft( Vector< complex<Type> > &xn );
        void ifft( Vector< complex<Type> > &Xk );
		void fft( const Vector<Type> &xn, Vector< complex<Type> > &Xk );
        void ifft( const Vector< complex<Type> > &Xk, Vector<Type> &xn );

	private:

        void bitReverse( Vector<int> &bitrev );
		void radix2( int nthpo, Complex<Type> *c0, Complex<Type> *c1 );
		void radix4( int nthpo, Complex<Type> *c0, Complex<Type> *c1,
                           Complex<Type> *c2, Complex<Type> *c3 );
		void radix8( int nxtlt, int nthpo, int length,
                     Complex<Type> *cc0, Complex<Type> *cc1,
                     Complex<Type> *cc2, Complex<Type> *cc3,
                     Complex<Type> *cc4, Complex<Type> *cc5,
                     Complex<Type> *cc6, Complex<Type> *cc7 );
		void dft( int direction, int n, Complex<Type> *b );
	};
	//	class FFTMR


	#include <fftmr-impl.h>

}
// namespace splab


#endif
// FFTMR_H
