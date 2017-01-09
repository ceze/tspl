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
 *                                fftpf.h
 *
 * Fast Fourier Transform with prime factor algorithm
 *
 * This class is designed for calculating discrete Fourier transform and
 * inverse discrete Fourier transform of 1D signals by using prime factor
 * algorithms. The signals can be arbitrary length. The largest prime factor
 * of n must be less than or equal to the constant PRIMEFACTOR defined below.
 *
 * The general idea is to factor the length of the DFT "n" into factors that are
 * efficiently handled by the routines. A number of short DFT's are implemented
 * with a minimum of arithmetical operations and using(almost) straight line
 * code resultingin very fast execution when the factors of n belong to this
 * set. Especially radix-10 is optimized. Prime factors, that are not in the
 * set of short DFT's are handled with direct evaluation of the DFP expression.
 *
 * The algorithm is modified from "xFFT.h" of "Pratical Fourier Transform
 * and C++ Implementation" written by Hequan Sun.
 *
 * Zhang Ming, 2010-09, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef FFTPF_H
#define FFTPF_H


#include <vector.h>

#define     PRIMEFACTOR     37
#define     PRIMEFACTORHALF (PRIMEFACTOR+1)/2
#define     PRIMECOUNT      20


namespace splab
{

    template<class Type>
    class FFTPF
    {

    public:

        FFTPF();
        ~FFTPF();

        void    fft( const Vector<Type> &xn, Vector< complex<Type> > &Xk );
        void    ifft( const Vector< complex<Type> > &Xk, Vector<Type> &xn );

        void    fft( const Vector< complex<Type> > &xn,
                     Vector< complex<Type> > &Xk );
        void    ifft( const Vector< complex<Type> > &Xk,
                      Vector< complex<Type> > &xn );

    private:

        bool    bAlloc;
        int     mNOldSize, mNFactor, nNewFactorSize, nHalfFactorSize,
                groupOffset, dataOffset, blockOffset, adr, groupNo,
                dataNo, blockNo, twNo;
        int     mSofarRadix[PRIMECOUNT], mActualRadix[PRIMECOUNT],
                mRemainRadix[PRIMECOUNT];

        Type    tPI, t2PI, tIRT2, tIRT3, omega, twRe, twIim;
        Type    *pftwRe, *pftgRe, *pfzRe, *pfvRe, *pfwRe, *pftwIm, *pftgIm,
                *pfzIm, *pfvIm, *pfwIm;
        Type    twiddleRe[PRIMEFACTOR], trigRe[PRIMEFACTOR],  zRe[PRIMEFACTOR],
                twiddleIm[PRIMEFACTOR], trigIm[PRIMEFACTOR],  zIm[PRIMEFACTOR],
                vRe[PRIMEFACTORHALF],   wRe[PRIMEFACTORHALF],
                vIm[PRIMEFACTORHALF],   wIm[PRIMEFACTORHALF];

        Type    c3_1,
                c5_1,  c5_2,   c5_3,  c5_4,   c5_5,
                c7_1,  c7_2,   c7_3,  c7_4,   c7_5,   c7_6,
                c9_2,  c9_3,   c9_4,  c9_5,   c9_6,   c9_7,   c9_8,  c9_9,
                c11_1, c11_2,  c11_3, c11_4,  c11_5,  c11_6,  c11_7, c11_8,
                c11_9, c11_10, c13_1, c13_2,  c13_3,  c13_4,  c13_5, c13_6,
                c13_7, c13_8,  c13_9, c13_10, c13_11, c13_12, c16_2, c16_3,
                c16_4, c16_5;

        Type    ttmp,
                t1_re,  t1_im,  t2_re,  t2_im,  t3_re,  t3_im,  t4_re,  t4_im,
                t5_re,  t5_im,  t6_re,  t6_im,  t7_re,  t7_im,  t8_re,  t8_im,
                t9_re,  t9_im,  t10_re, t10_im, t11_re, t11_im, t12_re, t12_im,
                t13_re, t13_im, t14_re, t14_im, t15_re, t15_im, t16_re, t16_im,
                t17_re, t17_im, t18_re, t18_im, t19_re, t19_im, t20_re, t20_im,
                t21_re, t21_im, t22_re, t22_im,
                m1_re,  m1_im,  m2_re,  m2_im,  m3_re,  m3_im,  m4_re,  m4_im,
                m5_re,  m5_im,  m6_re,  m6_im,  m7_re,  m7_im,  m8_re,  m8_im,
                m9_re,  m9_im,  m10_re, m10_im, m11_re, m11_im, m12_re, m12_im;

        void    releaseMem();
        void    allocateMem();
        void    factorize( int n, int &nFact, int *fact);
        void    primeSetup( int nPoints );
        void    permute( const Vector<Type> &xn, Vector< complex<Type> > &yn );
        void    permute( const Vector< complex<Type> > &xn, Vector< complex<Type> > &yn,
                         bool bTrans=true );
        void    initTrig( int radix );
        void    radix2( Type *aRe, Type *aIm );
        void    radix3( Type *aRe, Type *aIm );
        void    radix4( Type *aRe, Type *aIm );
        void    radix5( Type *aRe, Type *aIm );
        void    radix7( Type *aRe, Type *aIm );
        void    radix8( Type *aRe, Type *aIm );
        void    radix9( Type *aRe, Type *aIm );
        void    radix10( Type *aRe, Type *aIm );
        void    radix11( Type *aRe, Type *aIm );
        void    radix13( Type *aRe, Type *aIm );
        void    radix16( Type *aRe, Type *aIm );
        void    radixOther( int radix );
        void    twiddleFFT( int sofarRadix, int radix, int remainRadix,
                            Vector< complex<Type> > &yn );

    };
    // class FFTPF


    #include <fftpf-impl.h>

}
// namespace splab


#endif
//FFTPF_H
